/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-core/details/matcl_options_details.h"
#include "matcl-core/options/matcl_options.h"
#include "matcl-core/IO/disp_data_provider.h"
#include "matcl-core/IO/base_io.h"
#include "matcl-core/options/options_disp.h"

#include <iomanip>
#include <algorithm>

namespace matcl { namespace details
{
//-----------------------------------------------------------------------------------
//                              option_impl
//-----------------------------------------------------------------------------------
std::string option_impl::pretty_type_name(const std::string& type_name)
{
    const char* str_name    = typeid(std::string).name();
    if (type_name == str_name)
        return "string";
    else
        return type_name;
};

//-----------------------------------------------------------------------------------
//                              options_impl
//-----------------------------------------------------------------------------------
static options_impl::mutex_type g_mutex_global;

options_impl::mutex_type& options_impl::get_mutex_global()
{
    return g_mutex_global;
};

options_impl::options_impl()
{};
options_impl::~options_impl()
{};
options_impl::option_map& options_impl::get_options_predefined()
{
    static option_map m_options_predefined;
    return m_options_predefined;
};
void options_impl::error_unregistered_option(const std::string& opt_name)
{
    std::ostringstream os;    
    throw error::option_unregistered(opt_name);
}
void details::options_impl::register_option(const option& opt)
{
    std::unique_lock<mutex_type> lock(get_mutex_global());

    const std::string& option_name = opt.name();

    auto pos = get_options_predefined().find(option_name);

    if (pos == get_options_predefined().end())
        get_options_predefined().insert(pos, option_map::value_type(option_name, opt));
    else
        pos->second = opt;
};

void options_impl::set(const option& option_value)
{
    std::unique_lock<mutex_type> lock(m_mutex_local);

    const std::string& option_name = option_value.name();

    auto pos = m_options.find(option_name);

    if (pos == m_options.end())
        m_options.insert(pos, option_map::value_type(option_name, option_value));
    else
        pos->second = option_value;
};
void options_impl::remove(const option& option_value)
{
    std::unique_lock<mutex_type> lock(m_mutex_local);

    const std::string& option_name = option_value.name();

    auto pos = m_options.find(option_name);

    if (pos == m_options.end())
        return;
    else
        m_options.erase(pos);
};
void options_impl::set(const options_impl& other)
{
    for (auto pos = other.m_options.begin(); pos != other.m_options.end(); ++pos)
    {
        this->set(pos->second);
    };
};
void options_impl::remove(const options_impl& other)
{
    for (auto pos = other.m_options.begin(); pos != other.m_options.end(); ++pos)
    {
        this->remove(pos->second);
    };
};

class disp_map;

class option_disp : public option_visitor
{
    public:
        enum class value_type
        {
            value, default
        };

    private:
        value_type                  m_type;
        Integer                     m_width;
        align_type                  m_align;
        const disp_stream*          m_disp_stream;
        const disp_data_provider*   m_owner;
        std::string                 m_string;

    public:
        option_disp(const disp_data_provider* ow, const disp_stream* ds, Integer width, align_type at,
                    value_type vt)
            :m_width(width), m_disp_stream(ds), m_owner(ow), m_align(at), m_type(vt)
        {};

        std::string     string() const;

        virtual void    visit(const optional<Integer>& val, const option& opt) override;
        virtual void    visit(const optional<Real>& val, const option& opt) override;
        virtual void    visit(const optional<Complex>& val, const option& opt) override;
        virtual void    visit(const optional<std::string>& val, const option& opt) override;
        virtual void    visit(const optional<bool>& val, const option& opt) override;

        virtual void    visit(const optional<std::vector<Integer>>& val, const option& opt) override;
        virtual void    visit(const optional<std::vector<Real>>& val, const option& opt) override;
        virtual void    visit(const optional<std::vector<Complex>>& val, const option& opt) override;
        virtual void    visit(const optional<std::vector<std::string>>& val, const option& opt) override;
        virtual void    visit(const optional<std::vector<bool>>& val, const option& opt) override;

    private:
        void            set(const std::string& val);
        void            set_unregistered();

        template<class T>
        void            visit_type(const optional<T>& val, const option& opt);
        
        template<class T>
        std::string     to_string(const T& val);

        template<class T>
        std::string     to_string(const std::vector<T>& val);

        std::string     to_string(bool val);
};

class disp_map : public disp_data_provider
{
    private:
        using option_map    = std::map<std::string, option>;
        using iterator      = option_map::const_iterator;

    private:
        const option_map&   m_option_map;
        iterator            m_pos;
        iterator            m_pos_row_header;
        iterator            m_V1;
        iterator            m_V2;
        int                 m_col;
        int                 m_col_V1;
        int                 m_col_V2;
        bool                m_is_help;

    public:
        disp_map(const option_map& map, bool is_help)
            :m_option_map(map), m_is_help(is_help)
        {
        };

        disp_map(const disp_map&) = delete;
        disp_map& operator=(const disp_map&) = delete;

        virtual Integer     rows() const override   { return (Integer)m_option_map.size(); }
        virtual Integer     cols() const override   { return m_is_help? 2 : 3; };

        virtual void        begin() override;           //point to first element of the matrix
        virtual void        hold() override;            //store current point in variable V1
        virtual void        hold_column() override;     //store current point in variable V2
        virtual void        next_row() override;        //go to next row in given column
        virtual void        next_column() override;     //go to next column in given row
        virtual void        restore() override;         //go to the point stored in variable V1
        virtual void        restore_column() override;  //go to the point stored in variable V2

        virtual void        begin_row_headers() override;
        virtual void        next_row_header() override;

        virtual std::string get_value(const disp_stream* ds, Integer width, align_type at, 
                                      Integer r, Integer c) const override
        { 
            (void)r;
            iterator pos                = m_pos;
            const option& opt           = pos->second;

            if (c == 0)
            {
                option_disp vg(this, ds,width, at, option_disp::value_type::value);
                opt.accept(vg);

                return vg.string();
            }
            else
            {
                if (m_is_help == false)
                {
                    if (c == 1)
                    {
                        option_disp vg(this, ds,width, at, option_disp::value_type::default);
                        opt.accept(vg);

                        return vg.string();
                    }
                    else
                    {
                        std::string val     = opt.description();
                        return val;
                    };
                }
                else
                {
                    std::string val         = opt.description();
                    return val;
                };
            };
        };
        virtual std::string get_row_name(const disp_stream* orig, Integer r) const
        {
            (void)orig;
            (void)r;
            std::string lab = m_pos_row_header->first;
            return lab;
        };

        virtual std::string get_rows_label(const disp_stream* orig) const override
        {
            (void)orig;
            return "option";
        };
        
        virtual std::string get_col_name(const disp_stream* orig, Integer c) const
        {
            (void)orig;
            if (m_is_help == false)
            {
                if (c == 0)
                    return "value";
                else if (c == 1)
                    return "default";
                else
                    return "description";
            }
            else
            {
                if (c == 0)
                    return "default";
                else
                    return "description";
            };
        };
        virtual align_type get_align_row_header(const disp_stream* orig) const override
        {
            (void)orig;
            return align_type::left;
        };
        virtual align_type  get_align_col(const disp_stream* orig, Integer c) const override
        {
            (void)orig;
            if (m_is_help == false)
            {
                if (c == 0 || c == 1)
                    return align_type::right;
                else
                    return align_type::left;
            }
            else
            {
                if (c == 0)
                    return align_type::right;
                else
                    return align_type::left;
            };
        };
};

void disp_map::begin()
{
    m_pos   = m_option_map.begin();
    m_col   = 0;
};

void disp_map::hold()
{
    m_col_V1    = m_col;
    m_V1        = m_pos;
};

void disp_map::hold_column()
{
    m_col_V2    = m_col;
    m_V2        = m_pos;
};

void disp_map::next_row()
{
    ++m_pos;
};

void disp_map::next_column()
{
    ++m_col;
};

void disp_map::restore()
{
    m_col       = m_col_V1;
    m_pos       = m_V1;
};

void disp_map::restore_column()
{
    m_col       = m_col_V2;
    m_pos       = m_V2;
};

void disp_map::begin_row_headers()
{
    m_pos_row_header = m_option_map.begin();
}

void disp_map::next_row_header()
{
    ++m_pos_row_header;
};

std::string option_disp::string() const
{
    return m_string;
};

void option_disp::set(const std::string& val)
{
    m_string = val;
};

void option_disp::set_unregistered()
{
    m_string = "unregistered option";
};

template<class T>
void option_disp::visit_type(const optional<T>& val, const option& opt)
{
    T value;

    if (this->m_type == value_type::value)
    {
        if (!val)
        {
            set("?");
            return;
        };
    
        value = val.value();
    }
    else
    {
        try
        {
            value = options::get_default<T>(opt);
        }
        catch(...)
        {
            set_unregistered();
            return;
        };
    };

    std::string res = to_string(value);    
    set(res);
};

template<class T>
std::string option_disp::to_string(const T& val)
{
    std::string res = m_owner->to_string(m_disp_stream, m_width, val, m_align);
    return res;
};

template<class T>
std::string option_disp::to_string(const std::vector<T>& val)
{
    std::ostringstream os;
    os << "{";

    if (val.size() >= 1)
    {
        T elem  = val[0];
        os << m_owner->to_string(m_disp_stream, 0, elem, align_type::left);
    };

    for (size_t i = 1; i < val.size(); ++i)
    {
        T elem  = val[i];
        os << ", ";
        os << m_owner->to_string(m_disp_stream, 0, elem, align_type::left);
    };

    os << "}";

    std::string res = m_owner->to_string(m_disp_stream, m_width, os.str(), m_align);
    return res;
};

std::string option_disp::to_string(bool val)
{
    std::string res = (val == true)? "true" : "false";
    res             = m_owner->to_string(m_disp_stream, m_width, res, m_align);
    return res;
}

void option_disp::visit(const optional<Integer>& val, const option& opt)
{
    return visit_type<Integer>(val,opt);
};
void option_disp::visit(const optional<Real>& val, const option& opt)
{
    return visit_type<Real>(val,opt);
};
void option_disp::visit(const optional<Complex>& val, const option& opt)
{
    return visit_type<Complex>(val,opt);
};
void option_disp::visit(const optional<std::string>& val, const option& opt)
{
    return visit_type<std::string>(val,opt);
};
void option_disp::visit(const optional<bool>& val, const option& opt)
{
    return visit_type<bool>(val,opt);
};

void option_disp::visit(const optional<std::vector<Integer>>& val, const option& opt)
{
    return visit_type<std::vector<Integer>>(val,opt);
};
void option_disp::visit(const optional<std::vector<Real>>& val, const option& opt)
{
    return visit_type<std::vector<Real>>(val,opt);
};
void option_disp::visit(const optional<std::vector<Complex>>& val, const option& opt)
{
    return visit_type<std::vector<Complex>>(val,opt);
};
void option_disp::visit(const optional<std::vector<std::string>>& val, const option& opt)
{
    return visit_type<std::vector<std::string>>(val,opt);
};
void option_disp::visit(const optional<std::vector<bool>>& val, const option& opt)
{
    return visit_type<std::vector<bool>>(val,opt);
};

void options_impl::disp(const disp_stream_ptr& ds, const options& disp_options_0) const
{
    std::unique_lock<mutex_type> lock(m_mutex_local);

    disp_map dm(m_options, false);

    options disp_opt = disp_options_0;
    if (!disp_opt.get<bool>(matcl::opt::disp::display_zero()))
    {
        //force to display zero values
        disp_opt.set(matcl::opt::disp::display_zero(true));
    };

    return matcl::disp(dm, ds, disp_opt);
};

options_impl::ptr_type options_impl::copy() const
{
    ptr_type ret(new options_impl());
    ret->m_options = this->m_options;
    ret->m_notifier= this->m_notifier;
    return ret;
};

void options_impl::help(const disp_stream_ptr& ds, const options& disp_options_0)
{
    //this lock creates deadlock, if help is called after global
    //initialization, then everything should be ok.
    //std::unique_lock<mutex_type> lock(m_mutex_global);

    disp_map dm(get_options_predefined(), true);

    options disp_opt = disp_options_0;
    if (!disp_opt.get<bool>(matcl::opt::disp::display_zero()))
    {
        //force to display zero values
        disp_opt.set(matcl::opt::disp::display_zero(true));
    };

    return matcl::disp(dm, ds, disp_opt);
};

void options_impl::clear()
{
    std::unique_lock<mutex_type> lock(m_mutex_local);
    m_options.clear();
};

bool details::options_impl::has_option(const option& option_value) const
{
    std::unique_lock<mutex_type> lock(m_mutex_local);

    std::string name = option_value.name();

    if (m_notifier)
        m_notifier->report(name);

    auto pos = m_options.find(name);

    if (pos == m_options.end())
        return false;

    return pos->second.has_value();
}

void details::options_impl::install_notifier(const notifier_ptr& notif)
{
    std::unique_lock<mutex_type> lock(m_mutex_local);
    m_notifier = notif;
};

Integer options_impl::size() const
{
    std::unique_lock<mutex_type> lock(m_mutex_local);
    return (Integer)m_options.size();
};

//-----------------------------------------------------------------------------------
//                              option_impl
//-----------------------------------------------------------------------------------
class disp_map_option_help : public disp_data_provider
{
    private:
        const option_impl*  m_option;
        const option&       m_owner;

    public:
        disp_map_option_help(const option_impl* opt, const option& owner)
            :m_option(opt), m_owner(owner)
        {};

        virtual Integer     rows() const override   { return 0; }
        virtual Integer     cols() const override   { return 0; };

        virtual std::string get_value(const disp_stream*, Integer, align_type, 
                                      Integer, Integer) const override
        { 
            return "";
        };
        virtual void display_empty_matrix(const disp_stream*,line_printer&, Integer, Integer) const override
        {
            return;
        };
        virtual void end_display(const disp_stream*, line_printer&) const override
        {
            return;
        };

        virtual void start_display(const disp_stream* orig, line_printer& p) const override
        {
            std::string name        = m_option->name();
            std::string type        = m_option->get_option_type();
            std::string description = m_option->description();

            std::string m_value;
            std::string m_default;
            {
                option_disp vg(this, orig, 0, align_type::left, option_disp::value_type::value);
                m_option->accept(vg, m_owner);
                m_value             = vg.string();
            }
            {
                option_disp vg(this, orig, 0, align_type::left, option_disp::value_type::default);
                m_option->accept(vg, m_owner);
                m_default            = vg.string();
            }
            p.disp_empty_line();

            Integer header_w        = 10;
            {
                std::ostringstream os;
                os << std::setw(header_w) << "option : ";
                os << to_string(orig, 0, name, align_type::left);

                std::string str     = os.str();
                Integer term_w      = orig->get_terminal_width();
                Integer new_line_w  = orig->new_line_width(term_w);
                Integer line_w      = term_w - new_line_w;
                Integer str_w       = (Integer)str.size();
                line_w              = std::min(line_w,str_w);
                line_w              = std::max(line_w,0);

                p.disp_string(str, align_type::left);

                str                 = std::string(line_w, '-');

                p.disp_string(str, align_type::left);
            }
            {
                std::ostringstream os;
                os << std::setw(header_w) << "type : ";
                os << to_string(orig, 0, type, align_type::left);
                p.disp_string(os.str(), align_type::left);
            };
            {
                std::ostringstream os;
                os << std::setw(header_w) << "value : ";
                os << to_string(orig, 0, m_value, align_type::left);
                p.disp_string(os.str(), align_type::left);
            }
            {
                std::ostringstream os;
                os << std::setw(header_w) << "default : ";
                os << to_string(orig, 0, m_default, align_type::left);
                p.disp_string(os.str(), align_type::left);
            };
            {
                p.disp_empty_line();
                p.disp_string("description : ", align_type::left);
                p.disp_string(description, align_type::left);
            }
        };
};

void option_impl::help(const disp_stream_ptr& ds, const option& opt, const options& disp_options_0) const
{
    disp_map_option_help dm(this, opt);

    options disp_opt = disp_options_0;
    if (!disp_opt.get<bool>(matcl::opt::disp::display_zero()))
    {
        //force to display zero values
        disp_opt.set(matcl::opt::disp::display_zero(true));
    };

    return matcl::disp(dm, ds, disp_opt);
};

}};