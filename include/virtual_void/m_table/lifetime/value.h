#pragma once

#include <stdexcept>
#include <type_traits>

#include "../../erased/data/has_m_table/has_m_table.h"
#include "../../erased/lifetime/value_trait.h"

namespace virtual_void::m_table {
using value_data_ptr = erased::data::value_ptr<erased::data::with_m_table>;
}

namespace virtual_void::erased {
using namespace virtual_void;
template <>
struct data_trait<m_table::value_data_ptr> : value_trait<data::has_m_table> {};
}  // namespace virtual_void::erased

namespace virtual_void::m_table {

using value = erased::virtual_void<value_data_ptr>;
template <typename T>
using typed_value = erased::virtual_typed<T, value_data_ptr>;

static_assert(erased::is_virtual_void<value>);
static_assert(erased::is_virtual_void<typed_value<int>>);

}  // namespace virtual_void::m_table