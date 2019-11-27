#pragma once

#include <string>
#include <stdexcept>

void report_error_impl(const char* msg)
{
    throw std::runtime_error(msg? std::string(msg) : "");
};

extern "C"
{

void matcl_report_error(const char* msg)
{
    return report_error_impl(msg);    
};

};