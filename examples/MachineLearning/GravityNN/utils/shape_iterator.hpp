#pragma once

#include <vector>
#include <iterator>

class ShapeIter {
    class ShapeIterImpl;

public:
    ShapeIter(const std::vector<size_t>& shape) : shape(shape) {}

    ShapeIterImpl begin() const {
        return ShapeIterImpl(shape, true);
    }

    ShapeIterImpl end() const {
        return ShapeIterImpl(shape, false);
    }

private:
    const std::vector<size_t> shape;

    class ShapeIterImpl : public std::iterator<std::input_iterator_tag, std::vector<size_t>> {
    public:
        ShapeIterImpl(const std::vector<size_t>& shape, bool start): shape(shape) {
            if (start) {
                this->cur_vals = std::vector<size_t>(shape.size(), 0);
            } else {
                this->cur_vals = shape;
            }
        }

        ShapeIterImpl& operator++() {
            for (int i = cur_vals.size() - 1; i >= 0; --i) {
                if (++cur_vals[i] >= shape[i]) {
                    if (i == 0) {
                        cur_vals = shape;
                        break;
                    }
                    cur_vals[i] = 0;
                } else {
                    break;
                }
            }
            return *this;
        }

        bool operator!=(const ShapeIterImpl& other) const { 
            return (cur_vals != other.cur_vals) || (shape != other.shape);
        }

        const std::vector<size_t>& operator*() const {
            return cur_vals;
        }

    private:
        const std::vector<size_t> shape;
        std::vector<size_t> cur_vals;
    };
};