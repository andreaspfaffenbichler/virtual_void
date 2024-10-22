#pragma once

#include "holder.h"

namespace virtual_void::erased::data {

template <typename T, typename... ARGS>
auto make_shared_const(ARGS&&... args) {
  return std::make_shared<T const>(std::in_place, std::forward<ARGS>(args)...);
}

}  // namespace virtual_void::erased::data