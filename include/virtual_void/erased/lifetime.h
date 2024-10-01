#pragma once


//#include "../perfect_typeid_hash/index_table.h"
//
//namespace virtual_void
//{
////++utillities
//	template<typename ... Ts> struct overload : Ts ... { using Ts::operator() ...; };
//	template<class... Ts> overload(Ts...) -> overload<Ts...>;
//
//	template< typename... ARGS >
//	struct type_list
//	{
//		// see https://vittorioromeo.info/index/blog/cpp20_lambdas_compiletime_for.html
//		static void for_each( auto&& call )
//		{
//			( call.template operator()< ARGS >(), ... );
//		}
//	};	
//
//	template<int N, typename... Ts> using nth_type =
//		typename std::tuple_element<N, std::tuple<Ts...>>::type;
//
//	template< typename... Ts> using first = nth_type< 0, Ts...>;
////---utillities
//
////+++Forward
//class m_table;
//using m_table_t = m_table;
//template< typename CLASS > constexpr m_table* m_table_of();
//
//using m_table_const_void = std::pair< const m_table*, const void* >;
//using virtual_void = std::pair< const m_table*, void* >;
//
//template< typename P > auto to_virtual_void( const P* p ){ return m_table_const_void{ m_table_of< P >(), p }; }
//template< typename P > auto to_virtual_void( P* p ){ return m_table_const_void{ m_table_of< P >(), p }; }
//
//using typeid_const_void = std::pair< const std::type_info&, const void* >;
//using typeid_void = std::pair< const std::type_info&, void* >;
//
//template< typename P > auto to_typeid_void( const P* p ){ return typeid_const_void{ typeid( *p ), p }; }
//template< typename P > auto to_typeid_void( P* p ){ return typeid_void{ typeid( *p ), p }; }
//
//template< typename >  struct self_pointer;
//template<>  struct self_pointer< void * >		{ template< typename CLASS > using type = CLASS*; };
//template<>  struct self_pointer< const void * >	{ template< typename CLASS > using type = const CLASS*; };
//
//class type_info_dispatch;
////---Forward
//
////+++concepts
//template< typename DISPATCH, typename VOID >
//concept MtableDispatchableVoid = requires( const DISPATCH& void_ )
//{
//    { void_.data() }	-> std::convertible_to< VOID >;
//    { void_.m_table() } -> std::convertible_to< const m_table* >;
//};
////---concepts
//
////+++meta data and algorithms
//namespace class_hierarchy
//{
//	template< typename CLASS > struct class_;
//	
//	struct base { using bases_ = type_list<>; };
//	template< typename... BASES > struct bases{ using bases_ = type_list< BASES... >; };
//
//	template< typename CLASS, bool deep = true >
//	constexpr void visit_class( auto visitor )
//	{
//		visitor.template operator()< CLASS >();
//		using bases = class_< CLASS >::bases_;
//		bases::for_each( [ & ]< typename BASE >()
//		{ 
//			visitor.template operator()< CLASS, BASE >();
//			if constexpr( deep )
//			{
//				visit_class< BASE, deep >( visitor );
//			}
//		});
//	}
//	template< typename CLASSES, bool deep = true >
//	constexpr void visit_classes( auto visitor )
//	{
//		CLASSES::for_each( [ & ]< typename CLASS >()
//		{ 
//			visit_class< CLASS, deep >( visitor );
//		});
//	}
//
//	using bases_t =	std::vector< const std::type_info* >;
//	struct class_with_bases
//	{
//		const std::type_info* self;
//		m_table* m_table;
//		bases_t bases;
//	};
//	using classes_with_bases = std::map< std::type_index, class_with_bases >;
//
//	constexpr auto declare_visitor( classes_with_bases& registry )
//	{
//		return overload
//			{ [&]< typename C >				
//				{ 
//					registry[ typeid( C ) ].self = &typeid( C );  
//					registry[ typeid( C ) ].m_table = m_table_of< C >();  
//				}
//			, [&]< typename C, typename B >	{ registry[ typeid( C ) ].bases.emplace_back( &typeid( B ) ); }
//			};	
//	}
//	template< typename CLASS, bool deep = true >
//	constexpr nullptr_t declare( classes_with_bases& registry )
//	{
//		class_hierarchy::visit_class< CLASS, deep >( declare_visitor( registry ) );
//		return {};
//	}
//	template< typename CLASSES, bool deep = true >
//	constexpr nullptr_t declare_all( classes_with_bases& registry )
//	{
//		class_hierarchy::visit_classes< CLASSES, deep >( declare_visitor( registry ) );
//		return {};
//	}
//	template< typename CLASS >
//	constexpr auto declare_deep( classes_with_bases& registry )
//	{
//		return declare< CLASS, true >( registry );
//	}
//	template< typename CLASS >
//	constexpr auto declare_shallow( classes_with_bases& registry )
//	{
//		return declare< CLASS, false >( registry );
//	}
//
//	constexpr void visit_hierarchy( const std::type_index& self, const classes_with_bases& classes_with_bases, auto visitor );
//	inline constexpr void visit_bases( const bases_t& bases, const classes_with_bases& classes_with_bases, auto visitor )
//	{
//		for( auto base : bases )
//			visit_hierarchy( *base, classes_with_bases, visitor );
//	}
//	inline constexpr void visit_hierarchy( const std::type_index& self, const classes_with_bases& classes_with_bases, auto visitor )
//	{
//		auto found = classes_with_bases.find( self );
//		if( found == classes_with_bases.end() )
//			throw std::runtime_error("class not registered.");
//		visitor( *found->second.self );
//		visit_bases( found->second.bases, classes_with_bases, visitor );
//	}
//}
////---meta data and algorithms
//
////+++open method dispatch
//struct domain
//{
//	class_hierarchy::classes_with_bases		classes;
//	std::vector< type_info_dispatch* >		method_dispatches;	
//};
//
//class error;
//
//class m_table
//{
//public:
//	constexpr m_table( const std::type_info& type_info )
//		: type_info_( type_info )
//	{}
//	constexpr const std::type_info& type() const { return type_info_; }
//	using m_table_target_t = void(*)();
//	constexpr void set_method( int method_index, m_table_target_t target )
//	{
//		ensure_size( method_index + 1  );
//		table_[ method_index ] = target;
//	}
//	template< typename TRAGET >
//	constexpr void set_method( int method_index, TRAGET target )
//	{
//		auto m_table_target = reinterpret_cast< m_table_target_t >( target );
//		set_method( method_index, m_table_target );
//	}
//	constexpr void clear()
//	{
//		table_.clear();
//	}
//	constexpr auto at( int method_index ) const
//	{
//		return table_.at( method_index );
//	}
//	constexpr auto operator[]( int method_index ) const
//	{
//		return table_[ method_index ];
//	}
//private:
//	const std::type_info& type_info_;
//	std::vector< m_table_target_t > table_;
//	constexpr void ensure_size( std::size_t s )
//	{
//		if( table_.size() >= s )
//			return;
//		table_.insert( table_.end(), s  - table_.size(), nullptr );
//	}
//};
//template< typename CLASS > constexpr m_table* m_table_of()
//{
//	static m_table m_table_{ typeid( CLASS ) };
//	return &m_table_;
//}
//
//class error : public std::runtime_error
//{
//	using std::runtime_error::runtime_error;
//};
//
//class type_info_dispatch
//{
//public:
//	using dispatch_target_t = void(*)();
//private:
//	using entry_t = std::pair< const std::type_info*, dispatch_target_t >;
//	using method_table_t = std::vector< entry_t >; // faster than map, slower than hash_map 
//	method_table_t dispatchTable_;
//	dispatch_target_t default_ = reinterpret_cast< dispatch_target_t >( &throw_not_found );
//	const int m_table_index_ = -1;
//	using dispatch_target_index_t  = perfect_typeid_hash::index_table< dispatch_target_t >;
//	std::unique_ptr< dispatch_target_index_t > dispatch_target_index_;
//public:
//	type_info_dispatch() = default;
//	type_info_dispatch( domain& domain )
//		: m_table_index_( (int)domain.method_dispatches.size() )
//	{ 
//		domain.method_dispatches.push_back( this ); 
//	}
//	int m_table_index() const { return m_table_index_; }
//	static void throw_not_found(){ throw error( "no target specified." ); }
//	void define_default( auto f )
//	{
//		default_ = reinterpret_cast< dispatch_target_t >( f );
//	}
//	auto get_default() const
//	{
//		return default_;
//	}
//	auto override_erased( const std::type_info& register_type_info, auto f )
//	{
//		auto t = reinterpret_cast< dispatch_target_t >( f );
//		if( is_defined( register_type_info ) )
//			throw error( "Method for type already registered." );
//		dispatchTable_.insert( lower_bound( &register_type_info ), entry_t{ &register_type_info, t } );
//		return definition{};
//	}
//	template< typename TARGET = dispatch_target_t >
//	TARGET is_defined( const std::type_info& type_info ) const
//	{
//		auto found = lower_bound( &type_info );
//		if( found != dispatchTable_.end() && found->first == &type_info )
//			return reinterpret_cast< TARGET >( found->second );
//		return nullptr;
//	}
//	template< typename TARGET = dispatch_target_t >
//	TARGET lookup( const std::type_info& type_info ) const
//	{
//		if( !dispatch_target_index_ )
//			throw "Not yet sealed.";
//		return reinterpret_cast< TARGET >( dispatch_target_index_->at( &type_info ) );
//	}
//	void seal()
//	{
//		dispatch_target_index_ = std::make_unique< dispatch_target_index_t >( dispatchTable_, get_default() );
//	}
//	struct definition{};
//private:
//	auto lower_bound( const std::type_info* i ) const
//	{
//		return std::lower_bound
//			( dispatchTable_.begin(), dispatchTable_.end(), entry_t{ i, nullptr }
//			, []( const entry_t& i, const entry_t& v ){ return i.first < v.first; }
//			);
//	}
//};
//
//template< typename R, typename... ARGS >
//class method;
//
//template< typename R, typename... ARGS >
//class method< R( ARGS... ) >
//{
//	static_assert 
//		(	std::same_as< first< ARGS... >,	void* > 
//		||	std::same_as< first< ARGS... >,	const void* > 
//		); 
//public:
//	using dispatch_t = typename first< ARGS... >;
//	template< typename CLASS > using class_param_t = self_pointer< dispatch_t >::template type< CLASS >;
//	using param_t = std::pair< const std::type_info&, dispatch_t >;
//	using virtual_void_t = std::pair< const m_table*, dispatch_t >;
//	using erased_function_t = R(*)( ARGS... );
//private:
//	type_info_dispatch methodTable_;
//public:
//	method( domain& domain )
//		: methodTable_( domain )
//	{}
//	void define_default( auto f )
//	{
//		methodTable_->define_default( f );
//	}
//	auto override_erased( const std::type_info& ti, erased_function_t f ) { return methodTable_.override_erased( ti, f ); }
//	template< typename CLASS, typename FUNCTION >
//	auto override_( FUNCTION f )
//	{
//		auto fp = ensure_function_ptr< CLASS, ARGS... >( f );
//		return methodTable_.override_erased( typeid( CLASS ), fp );
//	}
//	template< typename... OTHER_ARGS >
//	R operator()( const std::type_info& type_info, dispatch_t dispatched, OTHER_ARGS&&... args ) const
//	{
//		auto f = methodTable_.lookup< erased_function_t >( type_info );
//		return f( dispatched, std::forward< OTHER_ARGS >( args )... );
//	}
//	template< typename... OTHER_ARGS >
//	R operator()( const param_t& param, OTHER_ARGS&&... args ) const
//	{
//		return (*this)( param.first, param.second, std::forward< OTHER_ARGS >( args )... );
//	}
//	template< typename... OTHER_ARGS >
//	R operator()( const virtual_void_t& param, OTHER_ARGS&&... args ) const
//	{
//		const m_table& m_table = *param.first;
//		auto erased_function = reinterpret_cast< erased_function_t >( m_table[ methodTable_.m_table_index() ] );
//		return (erased_function)( param.second, std::forward< OTHER_ARGS >( args )... );
//	}
//	template< typename CLASS, typename... OTHER_ARGS >
//	R operator()( CLASS* param, OTHER_ARGS&&... args ) const
//	{
//		return (*this)( to_typeid_void( param ), std::forward< OTHER_ARGS >( args )... );
//	}
//	template< typename CLASS, typename... OTHER_ARGS >
//	R operator()( const std::shared_ptr< CLASS >& param, OTHER_ARGS&&... args ) const
//	{
//		return (*this)( param.get(), std::forward< OTHER_ARGS >( args )... );
//	}
//	template< typename POINTER, typename... OTHER_ARGS >
//	R operator()( const POINTER& pointer, OTHER_ARGS&&... args ) const
//		requires MtableDispatchableVoid< POINTER, dispatch_t >
//	{
//		virtual_void_t param{ pointer.m_table(), pointer.data() }; 
//		return (*this)( param, std::forward< OTHER_ARGS >( args )... );
//	}
//	template< typename POINTER, typename... OTHER_ARGS >
//	R call( const POINTER& pointer, OTHER_ARGS&&... args ) const
//		requires MtableDispatchableVoid< POINTER, dispatch_t >
//	{
//		virtual_void_t param{ pointer.m_table(), pointer.data() }; 
//		return (*this)( param, std::forward< OTHER_ARGS >( args )... );
//	}
//	erased_function_t is_defined( const std::type_info& type_info ) const
//	{
//		return reinterpret_cast< erased_function_t >( methodTable_.is_defined( type_info ) );
//	}
//	template< typename C > auto is_defined() const
//	{
//		return is_defined( typeid( C ) );
//	}
//	void seal()
//	{
//		return methodTable_.seal();
//	}
//private:
//	template< typename CLASS, typename DISPATCH, typename... OTHER_ARGS >
//	static auto ensure_function_ptr( auto functor ) // if functor is a templated operator() from a stateless function object, instantiate it now!;
//	{
//		using functor_t = decltype( functor );
//		if constexpr( std::is_pointer_v< functor_t > )
//		{
//			return functor;
//		}
//		else
//		{
//								
//			return +[]( class_param_t< CLASS > self, OTHER_ARGS&&... args )->R
//			{
//				return functor_t{}( self, std::forward< OTHER_ARGS >( args )... );
//			};
//		}
//	}
//};
//
//template< typename R, typename... ARGS >
//class factory;
//
//template< typename R, typename... ARGS >
//class factory< R( ARGS...) >
//{
//public:
//	using factory_function_t = R(*)( ARGS... );
//private:
//	type_info_dispatch methodTable_;
//public:
//	auto override_erased( const std::type_info& ti, factory_function_t f ) { return methodTable_.override_erased( ti, f ); }
//	template< typename CLASS, typename FACTORY >
//	auto override_( FACTORY f )
//	{
//		auto fp = ensure_factory_ptr< CLASS >( f );
//		return methodTable_.override_erased( typeid( CLASS ), reinterpret_cast< factory_function_t >( fp ) );
//	}
//	R operator()( const std::type_info& type_info, ARGS&&... args ) const
//	{
//		auto f = methodTable_.lookup< factory_function_t >( type_info );
//		return f( std::forward< ARGS >( args )... );
//	}
//	template< typename CLASS > R operator()( ARGS&&... args ) const // to simplify tests!
//	{
//		return (*this)( typeid( CLASS ), std::forward< ARGS >( args )... );
//	}
//	factory_function_t is_defined( const std::type_info& type_info ) const
//	{
//		return methodTable_.is_defined< factory_function_t >( type_info );
//	}
//	template< typename C > auto is_defined() const
//	{
//		return is_defined( typeid( C ) );
//	}
//	void seal() { methodTable_.seal(); }
//private:
//	template< typename CLASS >
//	static auto ensure_factory_ptr( auto functor ) // if functor is a templated operator() from a stateless function object, instantiate it now!;
//	{
//		using functor_t = decltype( functor );
//		if constexpr( std::is_pointer_v< functor_t > )
//		{
//			return functor;
//		}
//		else
//		{
//			return +[]( ARGS&&... args )->R
//			{
//				return functor_t{}. template operator()< CLASS >( std::forward< ARGS >( args )... );
//			};
//		}
//	}
//};
////---open method dispatch
//
////+++open method algorithms
//template< typename CLASS >
//void fill_with_overload( auto& method, const auto& wrapper )
//{
//	if( !method.is_defined< CLASS >() )
//		method.override_< CLASS >( wrapper );
//}
//template< typename TYPE_LIST >
//void fill_with_overloads( TYPE_LIST, auto& method, const auto& wrapper )
//{
//	class_hierarchy::visit_classes< TYPE_LIST >( 
//		overload
//		{ [&]< typename C >				{ fill_with_overload< C >( method, wrapper ); }
//		, [&]< typename C, typename B >	{}
//		});
//}
//template< typename... CLASSES >
//void fill_with_overloads( auto& method, const auto& wrapper )
//{
//	fill_with_overloads( type_list< CLASSES... >{}, method, wrapper );
//}
//
//inline constexpr auto find_declared_in_bases( const class_hierarchy::classes_with_bases& registry, const class_hierarchy::bases_t& bases, const type_info_dispatch& method )
//{
//	typename type_info_dispatch::dispatch_target_t found = nullptr;
//	visit_bases( bases, registry, [ & ]( const std::type_info& base )
//		{
//			if( !found )
//				found = method.is_defined( base );				
//		});
//	return found;
//}
//inline void interpolate( const class_hierarchy::classes_with_bases& classes, type_info_dispatch* method )
//{
//	for( const auto& [ self, class_with_bases ] : classes )
//		if( !method->is_defined( *class_with_bases.self ) )
//			if( auto found = find_declared_in_bases( classes, class_with_bases.bases, *method ) )
//				method->override_erased( *class_with_bases.self, found );
//}
//inline void interpolate( const domain& domain )
//{
//	for( const auto& method : domain.method_dispatches )
//		interpolate( domain.classes, method );
//}
//inline void set_m_table( const class_hierarchy::class_with_bases& class_, const type_info_dispatch& method )
//{
//	auto target = method.is_defined( *class_.self );
//	if( !target )
//		target = method.get_default();
//	class_.m_table->set_method( method.m_table_index(), target );
//}
//inline void seal( domain& domain )
//{
//	for( const auto& method : domain.method_dispatches )
//		method->seal();
//}
//inline void fix_m_tables( const class_hierarchy::classes_with_bases& classes, const type_info_dispatch& method )
//{
//	for( const auto& class_ : classes )
//		set_m_table( class_.second, method );
//}
//inline void fix_m_tables( const domain& domain )
//{
//	for( const auto& method : domain.method_dispatches )
//		fix_m_tables( domain.classes, *method );
//}
//template< typename CLASSES > auto declare_classes( CLASSES, domain& domain )
//{
//	return class_hierarchy::declare_all< CLASSES >( domain.classes );
//}
//template< typename... CLASSES > auto declare_classes( domain& domain )
//{
//	return declare_classes( type_list< CLASSES... >{}, domain );
//}
//inline void build_m_tables( domain& domain )
//{
//	if( domain.classes.empty() )
//		throw error( "no classes declared." );
//	interpolate( domain );
//	seal( domain );
//	fix_m_tables( domain );
//}
////---open method algorithms
//
////+++erased cast
//template< template< typename > typename CONST, typename FOUND, typename FROM > auto erased_cast_implementation_( auto* from, const std::type_info& to )
//{
//	typename CONST< void >::type found = nullptr;
//	class_hierarchy::visit_class< FROM >
//		( overload
//		{ [&]< typename C >				
//			{ 
//				if( !found&& typeid( FROM ) == to ) 
//					found = static_cast< CONST< void >::type >( from ); 
//			}
//		, [&]< typename, typename B >	
//			{	
//				if( !found ) 
//					found = erased_cast_implementation_< CONST, FOUND, B >( static_cast< CONST< B >::type >( from ), to ); }
//		});
//	return found;
//}
//struct erased_const_cast_implementation
//{
//	template< typename T > struct const_ {	using type = const T*; }; 
//	template< typename FROM > auto operator()( const FROM* from, const std::type_info& to )
//	{
//		return erased_cast_implementation_< const_, typeid_const_void, FROM >( from, to );
//	}
//};
//struct erased_cast_implementation
//{
//	template< typename T > struct non_const_ { using type = T; }; 
//	template< typename FROM > auto operator()( FROM* from, const std::type_info& to )
//	{
//		return erased_cast_implementation_< non_const_, typeid_void, FROM >( from, to );	
//	}
//};
//using erased_const_cast_method = method< const void*( const void*, const std::type_info& to ) >;
//using erased_cast_method = method< void*( const void*, const std::type_info& to ) >;
//void fill_const_cast_for( auto classes, erased_const_cast_method& method )
//{
//	fill_with_overloads( classes, method, erased_const_cast_implementation{} );
//}
//template< typename... CLASSES >
//void fill_const_cast_for( erased_const_cast_method& method )
//{
//	fill_const_cast_for( type_list< CLASSES... >{}, method );
//}
//void fill_cast_for( auto classes, erased_cast_method& method )
//{
//	fill_with_overloads( classes, method, erased_cast_implementation{} );
//}
//template< typename... CLASSES >
//void fill_cast_for( erased_cast_method& method )
//{
//	fill_cast_for( type_list< CLASSES... >{}, method );
//}
//template< typename TO >
//auto cast_to( const erased_const_cast_method& cast, const auto& from )
//{
//    if( auto void_ = cast( from, typeid( std::remove_const_t< TO > ) ) )
//        return static_cast< TO* >( void_ );
//	return static_cast< TO* >( nullptr );
//}    
////---erased cast
//
////+++lifetime 
//template< typename T > class typed_shared_const;
//
//struct abstract_data
//{
//	const m_table_t* m_table_ = nullptr;
//	void* data_ = nullptr;
//	abstract_data( const m_table* m_table, void* data )
//		: m_table_( m_table )
//		, data_( data )
//	{}
//};
//template< typename M >
//struct concrete_data : abstract_data
//{
//	M m_;
//	concrete_data( M&& m )
//		: abstract_data( m_table_of< M >(), std::addressof( m_ ) )
//		, m_( std::move( m ) )
//	{}
//	template< typename... ARGS >
//	concrete_data( std::in_place_t,  ARGS&&... args )
//		: abstract_data( m_table_of< M >(), std::addressof( m_ ) )
//		, m_( std::forward< ARGS >( args )... )
//	{}
//};
//
//using shared_abstract_data_ptr = std::shared_ptr< abstract_data >;
//class shared_const
//{ 
//protected:
//	shared_abstract_data_ptr ptr_;
//public:
//	shared_const( const shared_abstract_data_ptr& ptr )
//		: ptr_( ptr )
//	{}
//    const void* data() const { return ptr_->data_; }
//    const std::type_info& type() const { return ptr_->m_table_->type(); }
//	const m_table* m_table() const { return ptr_->m_table_; };
//};
//static_assert( MtableDispatchableVoid< const shared_const, const void* > );
//
//template< typename T >
//class typed_shared_const : public shared_const
//{
//private:
//   typed_shared_const( shared_const&& ptr ) noexcept
//        : shared_const( std::move( ptr ) )
//    {}
// public:
//	using wrapped_type = T;
//	using shared_const::shared_const;
//	typed_shared_const( T&& v ) noexcept
//		: shared_const( std::make_shared< concrete_data< T > >( std::forward< T >( v ) ) )
//	{}
//	typed_shared_const( const T& v ) noexcept
//		: shared_const( std::make_shared< concrete_data< T > >( T{ v } ) )
//	{}
//	template< typename... ARGS > 
//	typed_shared_const( std::in_place_t, ARGS&&... args ) noexcept
//		: shared_const( std::make_shared< concrete_data< T > >( std::in_place, std::forward< ARGS >( args )... ) )
//	{}
//    template< typename T, typename... ARGS > friend typed_shared_const< T > make_shared_const( ARGS&&... args );
//    template< typename T > friend typed_shared_const< T > as( shared_const source );
//	template< typename DERIVED >
//	typed_shared_const( const typed_shared_const< DERIVED >& rhs ) noexcept
//		requires ( std::derived_from< DERIVED, T > && !std::same_as< DERIVED, T > )
//		: shared_const( rhs )
//	{}
//	template< typename DERIVED >
//	typed_shared_const( typed_shared_const< DERIVED >&& rhs ) noexcept
//		requires ( std::derived_from< DERIVED, T > && !std::same_as< DERIVED, T > )
//		: shared_const( std::move( rhs ) )
//	{}
//	template< typename DERIVED >
//	typed_shared_const& operator=( const typed_shared_const< DERIVED >& rhs ) noexcept
//		requires ( std::derived_from< DERIVED, T > && !std::same_as< DERIVED, T > )
//	{
//		typed_shared_const clone{ rhs };
//		swap( *this, clone );
//		return *this;
//	}
//	template< typename DERIVED >
//	typed_shared_const& operator=( typed_shared_const< DERIVED >&& rhs ) noexcept
//		requires ( std::derived_from< DERIVED, T > && !std::same_as< DERIVED, T > )
//	{
//		(*this) = std::move( rhs );
//		return *this;
//	}
//    friend void swap( typed_shared_const& rhs, typed_shared_const& lhs ) noexcept
//    {
//        using std::swap;
//        swap( rhs.ptr_, lhs.ptr_ );
//    }
//    const T& operator*() const noexcept  { return *static_cast< const T* >( data() ); }
//    const T* operator->() const noexcept { return  static_cast< const T* >( data() ); }
//};
//static_assert( MtableDispatchableVoid< const typed_shared_const< nullptr_t >, const void* > );
//template< typename T, typename... ARGS > typed_shared_const< T > make_shared_const( ARGS&&... args )
//{
//	return { std::in_place, std::forward< ARGS >( args )... };
//}
//template< typename T > typed_shared_const< T > as( shared_const source )
//{
//    if( source.type() != typeid( T ) )
//        throw error( "source is: " + std::string( source.type().name() ) + "." );
//    return typed_shared_const< T >{ std::move( source ) };
//}
//
//using unique_abstract_data_ptr = std::unique_ptr< abstract_data >;
//class unique
//{
//	unique_abstract_data_ptr ptr_;
//protected:
//	unique( unique_abstract_data_ptr&& ptr )
//		: ptr_( std::move( ptr ) )
//	{}
//public:
//    void* data() const { return ptr_->data_; }
//    const std::type_info& type() const { return ptr_->m_table_->type(); }
//	const m_table* m_table() const { return ptr_->m_table_; };
//};
//static_assert( MtableDispatchableVoid< const unique, void* > );
//static_assert( MtableDispatchableVoid< const unique, const void* > );
//
//template< typename T >
//class typed_unique : public unique
//{
//    template< typename T > friend auto as( unique&& source );
//	typed_unique( unique&& u ) noexcept
//		: unique( std::move( u ) ) 
//	{}
//public:
//	using unique::unique;
//	typed_unique( T&& v ) noexcept
//		: unique( std::make_unique< concrete_data< T > >( std::move( v ) ) )
//	{}
//	template< typename... ARGS > 
//	typed_unique( std::in_place_t, ARGS&&... args ) noexcept
//		: unique( std::make_unique< concrete_data< T > >( std::in_place, std::forward< ARGS >( args )... ) )
//	{}
//    T& operator*() const { return  *static_cast< T* >( data() ); }
//    T* operator->() const { return  static_cast< T* >( data() ); }
//};
//static_assert( MtableDispatchableVoid< const typed_unique< nullptr_t >, void* > );
//static_assert( MtableDispatchableVoid< const typed_unique< nullptr_t >, const void* > );
//
//template< typename T > auto as( unique&& source )
//{
//    if( source.type() != typeid( T ) )
//        throw error( "source is: " + std::string( source.type().name() ) + "." );
//    return typed_unique< T >{ std::move( source ) };
//}
//template< typename T, typename... ARGS > typed_unique< T > make_unique( ARGS&&... args )
//{
//	return { std::in_place, std::forward< ARGS >( args )... };
//}
////---lifetime 
//
//}