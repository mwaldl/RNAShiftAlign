#ifndef RNASHIFTALIGN_SHIFTMATRIX_M_HH
#define RNASHIFTALIGN_SHIFTMATRIX_M_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/**
 * @file shiftmatrix_m.hh
 * @brief Define M matrix class for shift alignments (with templated element type)
 *
 * This file defines the 4D dynamic programming matrix M used in the
 * Sankoff-style bi-alignment algorithm. The matrix stores optimal scores
 * for aligning prefixes of two RNA sequences with a bounded shift between
 * sequence and structure alignments.
 */

#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <tuple>
#include <string>

namespace RNAShiftAlign {

    /**
     * @brief 4D matrix class for shift bi-alignments
     *
     * Provides access via operator()(y1, y2, y3, y4) where:
     * - (y1, y2) represent positions in the sequence alignment U
     * - (y3, y4) represent positions in the structure alignment V
     * - Enforces constraint: |y1-y3| ≤ maxshift and |y2-y4| ≤ maxshift
     *
     * Memory layout: linearized 4D array with shift offset encoding
     * Storage: O(adim * bdim * (2*maxshift+1)^2)
     */
    template <class T>
    class ShiftMatrixM {
    public:
        typedef T elem_t; //!< type of elements
        typedef typename std::vector<elem_t>::size_type
            size_type; //!< size type (from underlying vector)
        typedef std::tuple<size_type, size_type, size_type, size_type>
            size_tuple_type; //!< type for tuple of sizes

    protected:
        std::vector<elem_t> mat_; //!< vector storing the matrix entries
        size_type adim_;          //!< first dimension (length of sequence a + 1)
        size_type bdim_;          //!< second dimension (length of sequence b + 1)
        size_type a_offsets_dim_; //!< third dimension (2*maxshift + 1)
        size_type b_offsets_dim_; //!< fourth dimension (2*maxshift + 1)
        size_type maxshift_; //!< maximum distance between y1 and y3 and between
                             //!< y2 and y4

        /**
         * @brief Computes address/index in 1D vector from 4D matrix indices
         *
         * @param y1 first index (position in sequence a for U)
         * @param y2 second index (position in sequence b for U)
         * @param y3 third index (position in sequence a for V)
         * @param y4 fourth index (position in sequence b for V)
         *
         * @return index in vector
         *
         * @throws std::out_of_range if indices violate bounds or maxshift constraint
         *
         * @note This method is used for all internal access to the vector mat_
         * @note Index layout: (y1, y2, offset_a, offset_b) where
         *       offset_a = y3 - y1 + maxshift, offset_b = y4 - y2 + maxshift
         */
        size_type
        addr(size_type y1, size_type y2, size_type y3, size_type y4) const {
            // Validate indices
            if (y1 >= adim_) {
                throw std::out_of_range(
                    "Index y1 is out of range [0,adim). "
                    "(y1:" +
                    std::to_string(y1) + ", adim:" + std::to_string(adim_) +
                    ")");
            }
            if (y2 >= bdim_) {
                throw std::out_of_range(
                    "Index y2 is out of range [0,bdim). "
                    "(y2:" +
                    std::to_string(y2) + ", bdim:" + std::to_string(bdim_) +
                    ")");
            }
            if (y3 >= adim_) {
                throw std::out_of_range(
                    "Index y3 is out of range [0,adim). "
                    "(y3:" +
                    std::to_string(y3) + ", adim:" + std::to_string(adim_) +
                    ")");
            }
            if (y4 >= bdim_) {
                throw std::out_of_range(
                    "Index y4 is out of range [0,bdim). "
                    "(y4:" +
                    std::to_string(y4) + ", bdim:" + std::to_string(bdim_) +
                    ")");
            }

            // Validate maxshift constraint
            size_type shift_a = (y1 > y3) ? (y1 - y3) : (y3 - y1);
            if (shift_a > maxshift_) {
                throw std::out_of_range(
                    "Index y3 violates maxshift constraint [y1-maxshift, y1+maxshift]. "
                    "(y3:" +
                    std::to_string(y3) + ", y1:" + std::to_string(y1) +
                    ", maxshift:" + std::to_string(maxshift_) + ")");
            }

            size_type shift_b = (y2 > y4) ? (y2 - y4) : (y4 - y2);
            if (shift_b > maxshift_) {
                throw std::out_of_range(
                    "Index y4 violates maxshift constraint [y2-maxshift, y2+maxshift]. "
                    "(y4:" +
                    std::to_string(y4) + ", y2:" + std::to_string(y2) +
                    ", maxshift:" + std::to_string(maxshift_) + ")");
            }

            // Compute linearized index
            return y1 * bdim_ * a_offsets_dim_ * b_offsets_dim_ +
                y2 * a_offsets_dim_ * b_offsets_dim_ +
                (y3 - y1 + maxshift_) * b_offsets_dim_ +
                (y4 - y2 + maxshift_);
        }

    public:
        /**
         * @brief Empty constructor
         *
         * Creates an empty matrix with all dimensions set to 0
         */
        ShiftMatrixM()
            : mat_(),
              adim_(0),
              bdim_(0),
              a_offsets_dim_(0),
              b_offsets_dim_(0),
              maxshift_(0) {}

        /**
         * @brief Construct with dimensions, optionally initialize from array
         *
         * @param adim first dimension of matrix (length of sequence a + 1)
         * @param bdim second dimension of matrix (length of sequence b + 1)
         * @param maxshift maximum distance between y1 and y3 and between y2 and y4
         * @param from pointer to array of elements (optional)
         *
         * @note if from is provided, initializes matrix from array
         */
        ShiftMatrixM(size_type adim,
                     size_type bdim,
                     size_type maxshift,
                     const elem_t *from = nullptr)
            : mat_(adim * bdim * (maxshift * 2 + 1) * (maxshift * 2 + 1)),
              adim_(adim),
              bdim_(bdim),
              maxshift_(maxshift),
              a_offsets_dim_(maxshift * 2 + 1),
              b_offsets_dim_(maxshift * 2 + 1) {
            if (from != nullptr) {
                std::copy(from,
                          from + adim_ * bdim_ * a_offsets_dim_ * b_offsets_dim_,
                          mat_.begin());
            }
        }

        /**
         * @brief Access maxshift parameter
         *
         * @return maximum shift between structure and sequence alignment
         */
        size_type
        get_maxshift() const {
            return maxshift_;
        }

        /**
         * @brief Access matrix dimensions
         *
         * @return tuple of dimensions (adim, bdim, a_offsets_dim, b_offsets_dim)
         */
        size_tuple_type
        sizes() const {
            return size_tuple_type(adim_, bdim_, a_offsets_dim_,
                                   b_offsets_dim_);
        }

        /**
         * @brief Resize all four dimensions
         *
         * @param adim first dimension
         * @param bdim second dimension
         * @param maxshift maximum shift between structure and sequence alignment
         *
         * @post Matrix is resized, existing content may be lost
         */
        void
        resize(size_type adim, size_type bdim, size_type maxshift) {
            adim_ = adim;
            bdim_ = bdim;
            maxshift_ = maxshift;
            a_offsets_dim_ = 2 * maxshift + 1;
            b_offsets_dim_ = 2 * maxshift + 1;
            mat_.resize(adim_ * bdim_ * a_offsets_dim_ * b_offsets_dim_);
        }

        /**
         * @brief Read-only access to matrix element
         *
         * @param y1 position in sequence a (for U)
         * @param y2 position in sequence b (for U)
         * @param y3 position in sequence a (for V)
         * @param y4 position in sequence b (for V)
         *
         * @return const reference to entry (y1,y2,y3,y4)
         */
        const elem_t &
        operator()(size_type y1,
                   size_type y2,
                   size_type y3,
                   size_type y4) const {
            return mat_[addr(y1, y2, y3, y4)];
        }

        /**
         * @brief Read/write access to matrix element
         *
         * @param y1 position in sequence a (for U)
         * @param y2 position in sequence b (for U)
         * @param y3 position in sequence a (for V)
         * @param y4 position in sequence b (for V)
         *
         * @return reference to entry (y1,y2,y3,y4)
         */
        elem_t &
        operator()(size_type y1, size_type y2, size_type y3, size_type y4) {
            return mat_[addr(y1, y2, y3, y4)];
        }

        /**
         * @brief Read access to matrix element (alternative interface)
         *
         * @param y1 position in sequence a (for U)
         * @param y2 position in sequence b (for U)
         * @param y3 position in sequence a (for V)
         * @param y4 position in sequence b (for V)
         *
         * @return copy of entry (y1,y2,y3,y4)
         */
        const elem_t
        get(size_type y1, size_type y2, size_type y3, size_type y4) const {
            return mat_[addr(y1, y2, y3, y4)];
        }

        /**
         * @brief Write access to matrix element (alternative interface)
         *
         * @param y1 position in sequence a (for U)
         * @param y2 position in sequence b (for U)
         * @param y3 position in sequence a (for V)
         * @param y4 position in sequence b (for V)
         * @param x element value to write
         */
        void
        set(size_type y1,
            size_type y2,
            size_type y3,
            size_type y4,
            const elem_t &x) {
            size_type idx = addr(y1, y2, y3, y4);
            if (idx >= mat_.size()) {
                throw std::out_of_range(
                    "ShiftMatrixM::set() computed index out of range: idx=" + std::to_string(idx) +
                    " >= size=" + std::to_string(mat_.size()) +
                    " for coordinates (" + std::to_string(y1) + "," + std::to_string(y2) + "," +
                    std::to_string(y3) + "," + std::to_string(y4) + ")"
                );
            }
            mat_[idx] = x;
        }

        /**
         * @brief Fill the whole matrix with the given value
         *
         * @param val value assigned to each entry
         * @post all matrix entries are set to val
         */
        void
        fill(const elem_t &val) {
            std::fill(mat_.begin(), mat_.end(), val);
        }

        /**
         * @brief Clear the matrix
         * @post the matrix is resized to dimensions (0,0,0,0)
         * @note behaves like std::vector::clear()
         */
        void
        clear() {
            adim_ = 0;
            bdim_ = 0;
            maxshift_ = 0;
            a_offsets_dim_ = 0;
            b_offsets_dim_ = 0;
            mat_.clear();
        }

        /**
         * @brief Transform matrix in place by applying a function to each element
         *
         * @param f function object or lambda
         *
         * @post All matrix entries are changed from x to f(x)
         * @note applies f via in-place std::transform to all matrix entries
         */
        template <class UnaryOperator>
        void
        transform(UnaryOperator f) {
            std::transform(mat_.begin(), mat_.end(), mat_.begin(), f);
        }

        /**
         * @brief Print matrix contents for debugging
         *
         * Outputs all entries in the M matrix including infinity values.
         * This helps identify which entries were not filled during computation.
         *
         * @param os Output stream
         * @param is_neg_infty Predicate to check if value is negative infinity
         */
        template <class Predicate>
        void
        debug_print(std::ostream& os, Predicate is_neg_infty) const {
            os << "=== M Matrix (4D: y1, y2, y3, y4) ===\n";
            os << "Dimensions: adim=" << adim_ << ", bdim=" << bdim_ << "\n";
            os << "max_shifts=" << maxshift_ << "\n\n";

            int count = 0;
            int infty_count = 0;
            for (size_type y1 = 0; y1 < adim_; ++y1) {
                for (size_type y2 = 0; y2 < bdim_; ++y2) {
                    for (size_type y3 = 0; y3 < adim_; ++y3) {
                        for (size_type y4 = 0; y4 < bdim_; ++y4) {
                            // Check delta_max constraint
                            size_type shift_a = (y1 > y3) ? (y1 - y3) : (y3 - y1);
                            size_type shift_b = (y2 > y4) ? (y2 - y4) : (y4 - y2);
                            if (shift_a > maxshift_ || shift_b > maxshift_)
                                continue;

                            const elem_t& val = mat_[addr(y1, y2, y3, y4)];
                            if (is_neg_infty(val)) {
                                os << "M(" << y1 << "," << y2 << "," << y3 << "," << y4 << ") = -INF\n";
                                ++infty_count;
                            } else {
                                os << "M(" << y1 << "," << y2 << "," << y3 << "," << y4 << ") = " << val << "\n";
                            }
                            ++count;
                        }
                    }
                }
            }

            os << "\nTotal entries: " << count << " (finite: " << (count - infty_count) << ", -INF: " << infty_count << ")\n";
        }
    };

} // end namespace RNAShiftAlign

#endif // RNASHIFTALIGN_SHIFTMATRIX_M_HH
