#ifndef RNASHIFTALIGN_SHIFTMATRIX_D_HH
#define RNASHIFTALIGN_SHIFTMATRIX_D_HH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* @file Define classes for D matrix for shift alignments (with templated element type)
 *
 * This file contains two classes:
 * - ShiftOffsetMatrix: 4D matrix storing shifts (x1,x2,y1,y2) at arc endpoints
 * - ShiftMatrixD: 6D D matrix = 2D base pair matrix × 4D offset matrix
 *
 * ShiftMatrixD is the higher-level class that uses ShiftOffsetMatrix.
 */
#include <iostream>
#include <vector>
#include <assert.h>

#include <algorithm>

#include <LocARNA/matrix.hh>
#include <LocARNA/basepairs.hh>

namespace LocARNA {

    /*
      Define classe for the offsets in the sequence alignment for the for
      indices that define the two matched arcs.  Provides access via operator
      (int,int,int,int)
    */
    template <class T>
    class ShiftOffsetMatrix {
    public:
        typedef T elem_t; //!< type of elements
        typedef typename std::vector<elem_t>::size_type
            size_type; //!< size type (from underlying vector)
        typedef std::tuple<size_type, size_type, size_type, size_type>
            size_tuple_type; //!< type for tuple of sizes

    protected:
        std::vector<elem_t> mat_; //!< vector storing the matrix entries
        size_type dim_;           //!< size per dimension
        size_type maxshift_; //!< maximum shift between sequence and structure
                             //!< alignment

        /**
         * Computes address/index in 1D vector from 4 shifts
         *
         * @param x1 shift at beginning of matched arc in first sequence
         * @param x2 shift at beginning of matched arc in second sequence
         * @param y1 shift at end of matched arc in first sequence
         * @param y2 shift at end of matched arc in first sequence
         *
         * @return index in vector
         * @note this method is used for all internal access to the vector mat_
         */
        size_type
        addr(int x1, int x2, int y1, int y2) const {
            assert(x1 >= -(int)maxshift_ &&
                   x1 <= (int)maxshift_);
            assert(x2 >= -(int)maxshift_ && x2 <= (int)maxshift_);
            assert(y1 >= -(int)maxshift_ && y1 <= (int)maxshift_);
            assert(y2 >= -(int)maxshift_ && y2 <= (int)maxshift_);
            // Row-major 4D addressing: ((x1 * dim + x2) * dim + y1) * dim + y2
            return (((x1 + maxshift_) * dim_ + (x2 + maxshift_)) * dim_ + (y1 + maxshift_)) * dim_ + (y2 + maxshift_);
        }

    public:
        /**
         * Empty constructor
         *
         */
        ShiftOffsetMatrix() : mat_(), dim_(0), maxshift_(0) {}

        /**
         * Construct with dimensions, optionally initialize from array
         *
         * @param maxshift maximum distance between i and k and between j and l
         * @param from pointer to array of elements
         *
         * @note if from given and !=0 initialize from array from
         *
         */
        ShiftOffsetMatrix(size_type maxshift, const elem_t *from = nullptr)
            : mat_((maxshift * 2 + 1) * (maxshift * 2 + 1) * (maxshift * 2 + 1) * (maxshift * 2 + 1)),
              dim_(maxshift * 2 + 1),
              maxshift_(maxshift) {
            if (from == nullptr)
                return;

            std::copy(from, from + dim_ * dim_ * dim_ * dim_,
                      mat_.begin());
        }

        /**
         * Access maxshift
         *
         * @return maximum shift between structure and sequence alignment
         */
        size_type
        get_maxshift() const {
            return maxshift_;
        }

        /**
         * Access size
         *
         * @return size of matrix as tuple of dimensions
         */
        size_tuple_type
        sizes() const {
            return size_tuple_type(dim_, dim_, dim_, dim_);
        }

        /**
         * Resize all four dimensions
         *
         * @param maxshift maximum shift between structure and sequence
         * alignment
         */
        void
        resize(size_type maxshift) {
            maxshift_ = maxshift;
            dim_ = 2 * maxshift + 1;
            mat_.resize(dim_ * dim_ * dim_ * dim_);
        }

        /**
         * Read access to matrix element
         *
         * @param x1 shift at beginning of matched arc in first sequence
         * @param x2 shift at beginning of matched arc in second sequence
         * @param y1 shift at end of matched arc in first sequence
         * @param y2 shift at end of matched arc in second sequence
         *
         * @return entry (x1,x2,y1,y2)
         */
        const elem_t &
        operator()(int x1,
                   int x2,
                   int y1,
                   int y2) const {
            return mat_[addr(x1, x2, y1, y2)];
        }

        /**
         * Read/write access to matrix element
         *
         * @param x1 shift at beginning of matched arc in first sequence
         * @param x2 shift at beginning of matched arc in second sequence
         * @param y1 shift at end of matched arc in first sequence
         * @param y2 shift at end of matched arc in second sequence
         *
         * @return reference to entry (x1,x2,y1,y2)
         */
        elem_t &
        operator()(int x1, int x2, int y1, int y2) {
            return mat_[addr(x1, x2, y1, y2)];
        }

        /**
         * Read access to matrix element
         *
         * @param x1 shift at beginning of matched arc in first sequence
         * @param x2 shift at beginning of matched arc in second sequence
         * @param y1 shift at end of matched arc in first sequence
         * @param y2 shift at end of matched arc in second sequence
         *
         * @return entry (x1,x2,y1,y2)
         */
        const elem_t
        get(int x1, int x2, int y1, int y2) const {
            return mat_[addr(x1, x2, y1, y2)];
        }

        /**
         * Write access to matrix element
         *
         * @param x1 shift at beginning of matched arc in first sequence
         * @param x2 shift at beginning of matched arc in second sequence
         * @param y1 shift at end of matched arc in first sequence
         * @param y2 shift at end of matched arc in second sequence
         * @param x element value
         *
         */
        void
        set(int x1,
            int x2,
            int y1,
            int y2,
            const elem_t &x) {
            mat_[addr(x1, x2, y1, y2)] = x;
        }

        /**
         * \brief Fill the whole matrix with the given value
         *
         * @param val value assigned to each entry
         * @post all matrix entries are set to val
         */
        void
        fill(const elem_t &val) {
            for (size_type i = 0; i < dim_ * dim_ * dim_ * dim_; ++i)
                mat_[i] = val;
        }

        /**
         * Clear the matrix
         * @post the matrix is resized to dimensions (0,0,0,0)
         * @note behaves like std::vector::clear()
         */
        void
        clear() {
            mat_.resize(0); // todo: is this needed?
            dim_ = 0;       // todo: is this needed?
            maxshift_ = 0;  // todo: is this needed?
            mat_.clear();
        }

        /**
         * Transform matrix in place due to applying a given function to each
         * element
         *
         * @param f function object
         *
         * @post All matrix entries are changed from x to f(x)
         * @note applies f via in place std::transform to all matrix entries
         */
        template <class UnaryOperator>
        void
        transform(UnaryOperator f) {
            std::transform(mat_.begin(), mat_.end(), mat_.begin(), f);
        }

        /**
         * @brief Print matrix contents for debugging
         *
         * Outputs all entries including infinity values to help identify
         * which entries were not filled during computation.
         *
         * @param os Output stream
         * @param is_neg_infty Function to check if value is negative infinity
         */
        template <class Predicate>
        void
        debug_print(std::ostream& os, Predicate is_neg_infty) const {
            os << "  === Offset Matrix (4D: x1, x2, y1, y2) ===\n";
            os << "  max_shifts=" << maxshift_ << "\n";

            int count = 0;
            int infty_count = 0;
            for (int x1 = -static_cast<int>(maxshift_); x1 <= static_cast<int>(maxshift_); ++x1) {
                for (int x2 = -static_cast<int>(maxshift_); x2 <= static_cast<int>(maxshift_); ++x2) {
                    for (int y1 = -static_cast<int>(maxshift_); y1 <= static_cast<int>(maxshift_); ++y1) {
                        for (int y2 = -static_cast<int>(maxshift_); y2 <= static_cast<int>(maxshift_); ++y2) {
                            const elem_t& val = mat_[addr(x1, x2, y1, y2)];
                            if (is_neg_infty(val)) {
                                os << "    (" << x1 << "," << x2 << "," << y1 << "," << y2 << ") = -INF\n";
                                ++infty_count;
                            } else {
                                os << "    (" << x1 << "," << x2 << "," << y1 << "," << y2 << ") = " << val << "\n";
                            }
                            ++count;
                        }
                    }
                }
            }

            os << "  Entries: " << count << " (finite: " << (count - infty_count) << ", -INF: " << infty_count << ")\n";
        }
    };

    /*
      Define classe for 6D D matrix for shift alignment.
    */
    template <class T>
    class ShiftMatrixD {
    public:
        typedef T elem_t; //!< type of elements
        typedef typename std::vector<elem_t>::size_type
            size_type; //!< size type (from underlying vector)
        typedef std::tuple<size_type, size_type, size_type>
            size_tuple_type; //!< type for tuple of sizes
        typedef BasePairs__Arc Arc;


    protected:
        Matrix<ShiftOffsetMatrix<elem_t>> mat_; //!< matrix storing the matrix entries
        size_type a_bps_dim_;           //!< size first dimension
        size_type b_bps_dim_;           //!< size second dimension
        size_type offset_dim_;          //!< size per offset dimension
        size_type maxshift_; //!< maximum shift between sequence and structure

    public:
        /**
         * Empty constructor
         *
         */
        ShiftMatrixD()
            : mat_(),
              a_bps_dim_(0),
              b_bps_dim_(0),
              offset_dim_(0),
              maxshift_(0) {}

        /**
         * Construct with dimensions
         *
         * @param maxshift maximum shift between structure and sequence of same
         * RNA
         *
         */
        ShiftMatrixD(size_type adim, size_type bdim, size_type maxshift)
            : mat_(adim, bdim),
              a_bps_dim_(adim),
              b_bps_dim_(bdim),
              offset_dim_(maxshift * 2 + 1),
              maxshift_(maxshift) {}

        /**
         * Access maxshift
         *
         * @return maximum shift between structure and sequence alignment
         */
        size_type
        get_maxshift() const {
            return maxshift_;
        }

        /**
         * Access size
         *
         * @return size of matrix as tuple of dimensions (number of bps in RNA
         * a, number of bps in RNA b, number of offsets considered)
         */
        size_tuple_type
        sizes() const {
            return size_tuple_type(a_bps_dim_, b_bps_dim_, offset_dim_);
        }

        /**
         * Resize both dimensions
         *
         * @param a_bps_dim first dimension
         * @param b_bps_dim second dimension
         * alignment
         */
        void
        resize(size_type a_bps_dim, size_type b_bps_dim) {
            a_bps_dim_ = a_bps_dim;
            b_bps_dim_ = b_bps_dim;
            mat_.resize(a_bps_dim, b_bps_dim);
        }

        /**
         * resize offset matrix for matched enclosing accoring to maxshift
         *
         * @param a index of bp in RNA a
         * @param b index of bp in RNA b
         */
        void
        create_offsetmatrix(size_type a, size_type b) {
            mat_(a, b).resize(maxshift_);
        }

        /**
         * Read access to matrix element by matched bps
         *
         * @param a index of bp in RNA a
         * @param b index of bp in RNA b
         *
         * @return entry (a,b)
         */
        const elem_t &
        operator()(size_type a, size_type b) const {
            return mat_(a, b);
        }

        /**
         * Read/write access to matrix element
         *
         * @param a index of bp in RNA a
         * @param b index of bp in RNA b
         *
         * @return reference to entry (a,b)
         */
        elem_t &
        operator()(size_type a, size_type b) {
            return mat_(a, b);
        }

        /**
         * Read access to matrix element by index
         *
         * @param matched bp in RNA a as Arc
         * @param matched bp in RNA b as Arc
         * @param x1 index at beginning of matched arc in first sequence
         * @param x2 index at beginning of matched arc in second sequence
         * @param y1 index at end of matched arc in first sequence
         * @param y2 index at end of matched arc in first sequence
         *
         * @return entry (a,b,x1,x2,y1,y2)
         */
        const elem_t
        get(Arc a, Arc b, size_type x1, size_type x2, size_type y1, size_type y2) const {
            return mat_(a.idx(), b.idx())(x1-a.left(), x2-b.left(), y1-a.right(), y2-b.right());
        }

        /**
         * Write access to matrix element
         *
         * @param matched bp in RNA a as Arc
         * @param matched bp in RNA b as Arc
         * @param x1 index at beginning of matched arc in first sequence
         * @param x2 index at beginning of matched arc in second sequence
         * @param y1 index at end of matched arc in first sequence
         * @param y2 index at end of matched arc in first sequence
         * @param x element value
         *
         */
        void
        set(Arc a, Arc b, size_type x1, size_type x2, size_type y1, size_type y2,
            const elem_t &x) {
            mat_(a.idx(), b.idx())(x1-a.left(), x2-b.left(), y1-a.right(), y2-b.right()) = x;
        }

        /**
         * @brief Print matrix contents for debugging
         *
         * Outputs all entries in the D matrix including infinity values
         * to help identify which entries were not filled during computation.
         *
         * @param os Output stream
         * @param is_neg_infty Function to check if value is negative infinity
         */
        template <class Predicate>
        void
        debug_print(std::ostream& os, Predicate is_neg_infty) const {
            os << "=== D Matrix (6D: arcA_idx, arcB_idx, x1, x2, y1, y2) ===\n";
            os << "Dimensions: a_bps=" << a_bps_dim_ << ", b_bps=" << b_bps_dim_ << "\n";
            os << "max_shifts=" << maxshift_ << "\n\n";

            int arc_pair_count = 0;
            for (size_type a = 0; a < a_bps_dim_; ++a) {
                for (size_type b = 0; b < b_bps_dim_; ++b) {
                    const auto& offset_mat = mat_(a, b);
                    // Check if offset matrix is allocated (dim > 0)
                    if (std::get<0>(offset_mat.sizes()) == 0)
                        continue;

                    os << "Arc pair (" << a << "," << b << "):\n";
                    offset_mat.debug_print(os, is_neg_infty);
                    ++arc_pair_count;
                }
            }

            os << "\nTotal arc pairs with allocated offset matrices: " << arc_pair_count << "\n";
        }
    };

} // end namespace LocARNA

#endif // RNASHIFTALIGN_SHIFTMATRIX_D_HH
