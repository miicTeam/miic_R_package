// Simple linear allocator with fixed buffer size.
// Suitable for contiguous container (std::vector) or related objects of POD
// types, with size determined and fixed on construction.
// Based on:
// https://gist.github.com/stoyannk/8122bbd9871fa7b5fd10
// https://docs.microsoft.com/en-us/cpp/standard-library/allocators

#ifndef MIIC_LINEAR_ALLOCATOR_H
#define MIIC_LINEAR_ALLOCATOR_H

#include <cassert>
#include <cstdlib>  // std::malloc, std::free, std::size_t
#include <memory>   // std::unique_ptr

namespace miic {
namespace utility {

namespace detail {

using std::size_t;

class LinearAllocator {
 public:
  LinearAllocator(size_t size)
      : m_ptr_{static_cast<char*>(std::malloc(size))},
        m_capacity_{size},
        m_free_space_{size} {}

  ~LinearAllocator() { std::free(m_ptr_); }

  void* Allocate(size_t size, unsigned alignment /* power of 2 */) {
    //assert((alignment & (alignment - 1)) == 0);  // assert power of 2
    auto current_p = static_cast<void*>(m_ptr_ + (m_capacity_ - m_free_space_));
    auto return_p = std::align(alignment, size, current_p, m_free_space_);
    // no space
    if (!return_p) {
      assert(false && "Linear allocator full!");
      return nullptr;
    }
    m_free_space_ -= size;

    return return_p;
  }

  void Deallocate() {
    // do nothing
  }

  void Reset(size_t freeSpace) { m_free_space_ = freeSpace; }

  size_t CurrentFreeSpace() const { return m_free_space_; }

 private:
  char* m_ptr_;
  size_t m_capacity_;
  size_t m_free_space_;
};

extern thread_local std::unique_ptr<LinearAllocator> li_alloc_ptr;

template <typename T>
struct TempStdAllocator {
 public:
  typedef T value_type;

 public:
  TempStdAllocator() = default;
  template <typename U>
  constexpr TempStdAllocator(const TempStdAllocator<U>&) noexcept {}

  template <typename U>
  bool operator==(const TempStdAllocator<U>&) const noexcept {
    return true;
  }
  template <typename U>
  bool operator!=(const TempStdAllocator<U>&) const noexcept {
    return false;
  }

  T* allocate(const size_t count) const {
    return reinterpret_cast<T*>(li_alloc_ptr->Allocate(
        unsigned(count * sizeof(T)), sizeof(size_t) == 4 ? 8 : 16));
  }

  void deallocate(T* const p, size_t) const noexcept {
    // do nothing
  }
};

struct TempAllocatorScope {
 public:
  TempAllocatorScope() : m_Space(li_alloc_ptr->CurrentFreeSpace()) {}

  ~TempAllocatorScope() { li_alloc_ptr->Reset(m_Space); }

 private:
  size_t m_Space;
};

}  // namespace detail
using detail::li_alloc_ptr;
using detail::LinearAllocator;
using detail::TempAllocatorScope;
using detail::TempStdAllocator;
}  // namespace utility
}  // namespace miic

#endif
