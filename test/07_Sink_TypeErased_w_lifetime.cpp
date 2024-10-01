﻿#include <iostream>
#include <string>
#include <any>
#include <functional>
#include <map>
#include <vector>
#include <utility>
#include <typeindex>
#include <memory>

#include "include/catch.hpp"

#include "../../include/std26/proxy.h"
#include "../../include/virtual_void/virtual_void.h"

namespace
{

namespace DB
{
    using FactoryFunction = std::function< virtual_void::shared_const( const std::string& ) >;
    using SinkFunction = std::function< void( const virtual_void::shared_const& ) >;

    struct System
    {
        std::map< std::string, FactoryFunction > factories;
        void Query( std::string what, const SinkFunction& sink )
        {
            auto dataLines = std::vector< std::pair< std::string, std::string > >  
            {   { "i", "1" }
            ,   { "i", "4711" }
            ,   { "s", "hello" }
            ,   { "s", "world" }
            ,   { "d", "3.14" }
            ,   { "ss", "chief" }
            };
            for( auto [ type, data ] : dataLines )
                sink( factories[ type ]( data ) );
        }        
    };
}

namespace Application
{
    struct IntData
    {
        int data = 0;
    }; 
    struct StringData
    {
        StringData( const std::string& s1 ) 
            : data( s1 )
        {}
        std::string data;
    }; 
    struct SuperStringData : StringData
    {
        SuperStringData( const std::string& s1, const std::string& s2 ) 
            : StringData( s1 )
            , more( s2 )
        {}
        std::string more;
    }; 
    struct DoubleData
    {
        double data;
    }; 

    std::string ToString_( const auto* x )
    {
        return std::to_string( x->data );
    }

    std::string ToString_( const StringData* x )
    {
        return x->data;
    }

    std::string ToString_( const SuperStringData* x )
    {
        return x->data + " -> " + x->more;
    }

    virtual_void::domain applicationDomain;

    auto entityToOut = virtual_void::method< void( const void* ) >{ applicationDomain };
    auto toString = virtual_void::method< std::string( const void* ) >{ applicationDomain };
    auto erased_const_ = virtual_void::erased_const_cast_method{ applicationDomain };

    void IntToOut( const IntData* i ){ std::cout << "int: " << i->data << std::endl; }

    void AnywhereInTheApplication()
    {
        entityToOut.override_< IntData >( &IntToOut );
    }
}

}
namespace virtual_void::class_hierarchy
{
	using namespace Application;
	template<> struct class_< StringData >: base {};
	template<> struct class_< SuperStringData >: bases< StringData >{};
	template<> struct class_< IntData >: base {};
	template<> struct class_< DoubleData >: base {};
}

namespace
{

namespace Application
{
    using classes = virtual_void::type_list< SuperStringData, StringData, IntData, DoubleData >;
}


TEST_CASE( "07_Sink_TypeErased_w_lifetime" ) 
{
    using namespace Application;
    using namespace DB;

    declare_classes( classes{}, applicationDomain );
 
    AnywhereInTheApplication();

    DB::System db;

    entityToOut.override_< StringData >( []( const StringData* s ) 
        { 
            std::cout << "string: " << s->data << std::endl; 
        });
 
    virtual_void::fill_with_overloads( classes{}, toString, []( const auto* x ){ return ToString_( x ); } );
    virtual_void::fill_const_cast_for( classes{}, erased_const_ );

    build_m_tables( applicationDomain );

    db.factories[ "i" ] = []( const std::string& data ){  return virtual_void::make_shared_const< IntData >( std::atoi( data.c_str() ) ); };
    db.factories[ "s" ] = []( const std::string& data ){  return virtual_void::make_shared_const< StringData >( data ); };
    db.factories[ "ss" ] = []( const std::string& data ){  return virtual_void::make_shared_const< SuperStringData >( data, "boss" ); };
    db.factories[ "d" ] = []( const std::string& data ){  return virtual_void::make_shared_const< DoubleData >( std::atof( data.c_str() ) ); };

    db.Query( "junk", []( const virtual_void::shared_const& e )
    { 
        // call open method
        std::cout << "type_info: " << e.type().name() << ": " << toString( e ) << std::endl;
        
        try
        {
            // call open method
            entityToOut( e );
        }
        catch( std::exception& e )
        {
            std::cout << "error: " << e.what() << std::endl;
        }

        // cast back from erased -> "unerase"
        if( auto stringData = cast_to< const StringData >( erased_const_, e ) )
            std::cout << "stringData: " << stringData->data << std::endl;
    });

 
}

}