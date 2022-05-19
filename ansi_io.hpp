#ifndef _ANSI_IO_HPP_
#define _ANSI_IO_HPP_

namespace ansi
{
	const string red     = "\x1b[1;31m";
	const string green   = "\x1b[1;32m";
	const string yellow  = "\x1b[1;33m";
	const string blue    = "\x1b[1;34m";
	const string magenta = "\x1b[1;35m";
	const string cyan    = "\x1b[1;36m";
	const string white   = "\x1b[1;37m";
	const string reset   = "\x1b[0m";
	const string erase   = "\x1b[2K\x1b[0m";
}

#endif
