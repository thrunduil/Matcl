#pragma once

#include <streambuf>
#include <ios>
#include <algorithm>
#include <iterator>

namespace matcl { namespace details
{

template <class FwIter>
class rangebuf : public std::streambuf 
{
	private:
		FwIter		first;
		FwIter		last;
		FwIter		current;

		static const int bufsize = 128;

		char		buf[bufsize];

		void		fill_buffer();

	public:
		rangebuf(FwIter f, FwIter l);

	protected:
		pos_type	seekoff(off_type off, std::ios::seekdir way, std::ios::openmode which);
		pos_type	seekpos(pos_type pos, std::ios::openmode which);
  
		int			underflow();
		int			pbackfail(int c);
};

template <class FwIter>
rangebuf<FwIter>::rangebuf(FwIter f, FwIter l)
  : std::streambuf(), first(f), last(l), current(f)
{
	fill_buffer();
}

template <class FwIter>
void rangebuf<FwIter>::fill_buffer()
{
	ptrdiff_t len = std::distance(current, last);
	ptrdiff_t n = std::min(len, ptrdiff_t(bufsize));
	FwIter tmp = current;
	std::advance(tmp, n);
	std::copy(current, tmp, buf);
	setg(buf, buf, buf + n);
}

template <class FwIter>
int rangebuf<FwIter>::underflow() 
{
	std::advance(current, egptr() - eback());
	fill_buffer();
	return (gptr() != egptr()) ? static_cast<int>(static_cast<unsigned char>(*gptr())) : EOF;
}

template <class FwIter>
int rangebuf<FwIter>::pbackfail(int c) 
{
	if (gptr() == eback() && current != first) 
	{
		FwIter tmp = first;
		std::advance(tmp, std::distance(first, current) - 1);
		current = tmp;
		fill_buffer();
		gbump(1);
	}

	if (gptr() == eback() || (c != EOF && c != *(gptr() - 1)))
	{
		return EOF;
	}
	else 
	{
		gbump(-1);
		return static_cast<unsigned char>(*gptr());
	}
}

template <class FwIter>
std::streambuf::pos_type
rangebuf<FwIter>::seekoff(off_type off, std::ios::seekdir way,
                          std::ios::openmode which) 
{
	switch(way) 
	{
		case std::ios::beg:
			return seekpos(pos_type(off), which);
		case std::ios::cur:
			return seekpos(pos_type(off + std::distance(first, current)), which);
		case std::ios::end:
			return seekpos(pos_type(off + std::distance(first, last)), which);
		default:
			return pos_type(off_type(-1));
	}
}

template <class FwIter>
std::streambuf::pos_type
rangebuf<FwIter>::seekpos(pos_type pos, std::ios::openmode which) 
{
	if (which != std::ios::in)
	    return pos_type(off_type(-1));

	off_type offset = off_type(pos);

	if (offset < 0 || offset > std::distance(first, last))
		return pos_type(off_type(-1));

	FwIter tmp = first;
	std::advance(tmp, offset);
	current = tmp;
	fill_buffer();
	
	return pos;
}

}};