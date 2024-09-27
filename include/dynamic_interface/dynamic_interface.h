#include <memory>
#include <type_traits>

namespace dynamic_interface
{
    template< typename ERASED >
    struct trait; 

    template<>
    struct trait< void* >
    {
        using type = void*;

        using param_t = void*;

        template< typename FROM >
        static void* erase( FROM&& from )
        {
            return static_cast< std::remove_cvref_t< FROM > * >( &from );
        }
        template< typename TO >
        static auto unerase( void* from )
        {
            return static_cast< std::remove_cvref_t< TO > * >( from );
        }
    };
};

#define _detail_EXPAND(...) _detail_EXPAND4(_detail_EXPAND4(_detail_EXPAND4(_detail_EXPAND4(__VA_ARGS__))))
#define _detail_EXPAND4(...) _detail_EXPAND3(_detail_EXPAND3(_detail_EXPAND3(_detail_EXPAND3(__VA_ARGS__))))
#define _detail_EXPAND3(...) _detail_EXPAND2(_detail_EXPAND2(_detail_EXPAND2(_detail_EXPAND2(__VA_ARGS__))))
#define _detail_EXPAND2(...) _detail_EXPAND1(_detail_EXPAND1(_detail_EXPAND1(_detail_EXPAND1(__VA_ARGS__))))
#define _detail_EXPAND1(...) __VA_ARGS__

#define _detail_EXPAND_(...) _detail_EXPAND_4(_detail_EXPAND_4(_detail_EXPAND_4(_detail_EXPAND_4(__VA_ARGS__))))
#define _detail_EXPAND_4(...) _detail_EXPAND_3(_detail_EXPAND_3(_detail_EXPAND_3(_detail_EXPAND_3(__VA_ARGS__))))
#define _detail_EXPAND_3(...) _detail_EXPAND_2(_detail_EXPAND_2(_detail_EXPAND_2(_detail_EXPAND_2(__VA_ARGS__))))
#define _detail_EXPAND_2(...) _detail_EXPAND_1(_detail_EXPAND_1(_detail_EXPAND_1(_detail_EXPAND_1(__VA_ARGS__))))
#define _detail_EXPAND_1(...) __VA_ARGS__
#define _detail_RMCVREF(x) typename std::remove_const<typename std::remove_volatile<typename std::remove_reference<x>::type>::type>::type
#define _detail_PARENS ()
#define _detail_foreach_macro_h(macro, a, ...) macro(a) \
__VA_OPT__(_detail_foreach_macro_a _detail_PARENS (macro, __VA_ARGS__))
#define _detail_foreach_macro_a() _detail_foreach_macro_h
#define _detail_foreach_macro(macro, ...) _detail_EXPAND(_detail_foreach_macro_h(macro, __VA_ARGS__))
#define _detail_map_macro_h(macro, a, ...) macro(a) \
__VA_OPT__(, _detail_map_macro_a _detail_PARENS (macro, __VA_ARGS__))
#define _detail_map_macro(macro, ...) _detail_EXPAND(_detail_map_macro_h(macro, __VA_ARGS__))
#define _detail_map_macro_a() _detail_map_macro_h
#define _detail_CONCAT_H(a, b) a ## b
#define _detail_CONCAT(a, b) _detail_CONCAT_H(a, b)
#define _detail_PARAM_LIST_H(b, c, f, ...) std::forward<decltype(c)>(c) __VA_OPT__(, _detail_PARAM_LIST_A _detail_PARENS (b, _detail_CONCAT(b, c), __VA_ARGS__))
#define _detail_PARAM_LIST_A() _detail_PARAM_LIST_H
#define _detail_PARAM_LIST(...) _detail_EXPAND_(_detail_PARAM_LIST_H(__VA_ARGS__))
#define _detail_PARAM_LIST_2H(b, c, f, ...) f c __VA_OPT__(, _detail_PARAM_LIST_2A _detail_PARENS (b, _detail_CONCAT(b, c), __VA_ARGS__))
#define _detail_PARAM_LIST_2A() _detail_PARAM_LIST_2H
#define _detail_PARAM_LIST2(...) _detail_EXPAND_(_detail_PARAM_LIST_2H(__VA_ARGS__))
#define _detail_EXPAND_LIST(...) __VA_ARGS__

#define _detail_INTERFACE_FUNCTION_PTR_DECL(type, name, ...) type (* name)(erased_param_t __VA_OPT__(, __VA_ARGS__));
#define _detail_LEAD_COMMA_H(...) __VA_OPT__(,)
#define _detail_INTERFACE_FPD_H(l) _detail_INTERFACE_FUNCTION_PTR_DECL l
#define _detail_INTERFACE_LIMP_H(l) _detail_INTERFACE_LAMBDA_IMPL l
#define _detail_INTERFACE_METHOD_H(l) _detail_INTERFACE_METHOD l
#define _detail_LEAD_COMMA_H_E(l) _detail_LEAD_COMMA_H l
#define _detail_INTERFACE_LAMBDA_IMPL(type, name, ...) \
name([](erased_param_t _vp __VA_OPT__(,_detail_PARAM_LIST2(a, _sig, __VA_ARGS__))) \
{return dynamic_interface::trait<erased_t>::unerase<_tp>(_vp)->name(__VA_OPT__(_detail_PARAM_LIST(a, _sig, __VA_ARGS__)));})
#define _detail_INTERFACE_METHOD(type, name, ...) \
type name(__VA_OPT__(_detail_PARAM_LIST2(a, _sig, __VA_ARGS__))) { \
    return _body.name(_body._ref __VA_OPT__(, _detail_PARAM_LIST(a, _sig, __VA_ARGS__))); \
}
#define _detail_DECLARE_INTERFACE( _erased, n, limp, l) \
class n { \
    using erased_t = _erased; \
    using erased_param_t = dynamic_interface::trait<_erased>::param_t; \
    struct _impl { \
        erased_t _ref = nullptr; \
        _detail_foreach_macro(_detail_INTERFACE_FPD_H, _detail_EXPAND_LIST l)\
        _impl() = default; \
        _impl(_impl&) = default;\
        _impl(_impl&&) = default;\
        template <typename _tp> \
        _impl(_tp&& v)  \
        : _ref(dynamic_interface::trait<_erased>::erase(std::forward<_tp>(v))) _detail_LEAD_COMMA_H_E(l) _detail_map_macro(limp, _detail_EXPAND_LIST l) {}\
    } _body;\
    public: \
    template <typename _tp> \
    n(_tp&& v) : _body(v) {} \
    _detail_foreach_macro(_detail_INTERFACE_METHOD_H, _detail_EXPAND_LIST l)    \
    n(const n&) = default;\
    n(n&) = default;\
    n(n&&) = default;\
};
#define DECLARE_INTERFACE(_erased, name, ...) _detail_DECLARE_INTERFACE(_erased, name, _detail_INTERFACE_LIMP_H, (__VA_ARGS__))
#define INTERFACE_METHOD(...) (__VA_ARGS__),


/*
THIS INTERFACE:
DECLARE_INTERFACE(example,
    (void, print, const char *)
)
EXPANDS TO:
    class example {
        struct _impl {
            void *_ref = nullptr;
            void (*print)(void *, const char *);
            _impl() = default;
            template <typename _tp>
            _impl(_tp &&v)
                : _ref(const_cast<
                    typename std ::remove_const<typename std ::remove_volatile<
                        typename std ::remove_reference<_tp>::type>::type>::type
                        *>(&v)),
                print([](void *_vp, const char *_sig) {
                    return static_cast<typename std ::remove_const<
                        typename std ::remove_volatile<
                            typename std ::remove_reference<_tp>::type>::type>::
                                            type *>(_vp)
                        ->print(std ::forward<decltype(_sig)>(_sig));
                }) {}
        } _body;

    public:
        example() = default;
        template <typename _tp> example(_tp &&v) : _body(v) {}
        void print(const char *_sig) {
            return _body.print(_body._ref, std ::forward<decltype(_sig)>(_sig));
        }
        operator bool() { return _body._ref != nullptr; }
    };
*/
