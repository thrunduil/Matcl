/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#pragma once

#include <map>
#include <string>
#include <iostream>

// A helper macro that wraps test::test_selector() method.
// Used to shorten calls to test_selector in lists of functions
//
// layer   The layer to provide to test_selector call.
// name    The name of test to provide to test_selector call.

#define SELECT_TEST(layer, name)                        \
    if (test_selector().is_selected((layer),#name))     \
           {   (name);     }                                   
                                                      

namespace matcl { namespace test
{

namespace details
{

    // Test selector implementation.
    class test_selector_impl
    {
        public:
            test_selector_impl()
            {};

            // Adds a selection of desired tests on layer.
            // Will add a test on a given layer. This limits testing to this particular case
            // on the layer
            //     layer   The layer to which to add.
            //     name    The name of test to restrict to.
            void add_selection(int layer, const std::string& name)
            {
                m_selection_map[layer] = name;
            };

            // Query if test is selected on layer
            //     layer   The layer to which to add.
            //     name    The name of test to restrict to.
            // return true if selected or layer was not initialized. False otherwise
            bool is_selected(int layer, const std::string& name)
            {
                return m_selection_map.count(layer) == 0 || 
                    m_selection_map[layer] == name;
            };

        private:
            std::map<int, std::string> m_selection_map;
    };
};

// Interface to static test selector.
class test_selector
{
    public:
        // Interface to add a selection of desired tests on layer.
        // Will add a test on a given layer. This limits testing to this 
        // particular case on the layer
        //     layer   The layer to which to add.
        //     name    The name of test to restrict to.
        void add_selection(int layer, const std::string& name)
        {
            selector.add_selection(layer, name);
        };

        // Interface to query if test is selected on layer
        //     layer   The layer to which to add.
        //     name    The name of test to restrict to.
        //
        // return true if selected or layer was not initialized. False otherwise
        bool is_selected(int layer, const std::string& name)
        {
            bool val = selector.is_selected(layer, name);
            return val;
        };
        
        // Just another version of is_selected invisible to scons' 
        // get_test_dictionary_from_sources
        bool is_selected_no_scons_parse(int layer, const std::string& name)
        {
            return is_selected(layer, name);
        };

    private:
        static details::test_selector_impl selector;
};

}};