#include "linear_allocator.h"

namespace miic {
namespace utility {
namespace detail {
thread_local std::unique_ptr<LinearAllocator> li_alloc_ptr;
}
}  // namespace utility
}  // namespace miic
