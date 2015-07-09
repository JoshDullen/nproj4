namespace Free.Ports.Proj4.LibcStuff
{
	public enum ERRORNUMBER
	{
		EGENERIC=0,
		EPERM,
		ENOENT,
		ESRCH,
		EINTR,
		EIO,
		ENXIO,
		E2BIG,
		ENOEXEC,
		EBADF,
		ECHILD,
		EAGAIN,
		ENOMEM,
		EACCES,
		EFAULT,
		ENOTBLK,
		EBUSY,
		EEXIST,
		EXDEV,
		ENODEV,
		ENOTDIR,
		EISDIR,
		EINVAL,
		ENFILE,
		EMFILE,
		ENOTTY,
		ETXTBSY,
		EFBIG,
		ENOSPC,
		ESPIPE,
		EROFS,
		EMLINK,
		EPIPE,
		EDOM,
		ERANGE,
		EDEADLK,
		EDEADLOCK,
		ERROR37,
		ENAMETOOLONG,
		ENOLCK,
		ENOSYS,
		ENOTEMPTY,
		EILSEQ,
		ERROR43,
		ERROR44,
		ERROR45,
		ERROR46,
		ERROR47,
		ERROR48,
		ERROR49,
		EPACKSIZE,
		EOUTOFBUFS,
		EBADIOCTL,
		EBADMODE,
		EWOULDBLOCK,
		EBADDEST,
		EDSTNOTRCH,
		EISCONN,
		EADDRINUSE,
		ECONNREFUSED,
		ECONNRESET,
		ETIMEDOUT,
		EURG,
		ENOURG,
		ENOTCONN,
		ESHUTDOWN,
		ENOCONN,
		EAFNOSUPPORT,
		EPROTONOSUPPORT,
		EPROTOTYPE,
		EINPROGRESS,
		EADDRNOTAVAIL,
		EALREADY,
		EMSGSIZE,
		ENOTSOCK,
		ENOPROTOOPT,
		ERROR76,
		ERROR77,
		ERROR78,
		ERROR79,
		STRUNCATE
	}

	public static partial class Libc
	{
		static readonly string[] _sys_errlist=new string[]
		{
			"Error 0",							// EGENERIC
			"Not owner",						// EPERM
			"No such file or directory",		// ENOENT
			"No such process",					// ESRCH
			"Interrupted system call",			// EINTR
			"I/O error",						// EIO
			"No such device or address",		// ENXIO
			"Arg list too long",				// E2BIG
			"Exec format error",				// ENOEXEC
			"Bad file number",					// EBADF
			"No children",						// ECHILD
			"Resource temporarily unavailable",	// EAGAIN
			"Not enough core",					// ENOMEM
			"Permission denied",				// EACCES
			"Bad address",						// EFAULT
			"Block device required",			// ENOTBLK
			"Resource busy",					// EBUSY
			"File exists",						// EEXIST
			"Cross-device link",				// EXDEV
			"No such device",					// ENODEV
			"Not a directory",					// ENOTDIR
			"Is a directory",					// EISDIR
			"Invalid argument",					// EINVAL
			"File table overflow",				// ENFILE
			"Too many open files",				// EMFILE
			"Not a typewriter",					// ENOTTY
			"Text file busy",					// ETXTBSY
			"File too large",					// EFBIG
			"No space left on device",			// ENOSPC
			"Illegal seek",						// ESPIPE
			"Read-only file system",			// EROFS
			"Too many links",					// EMLINK
			"Broken pipe",						// EPIPE
			"Math argument",					// EDOM
			"Result too large",					// ERANGE
			"Resource deadlock avoided",		// EDEADLK
			"Resource deadlock avoided",		// EDEADLOCK
			"Unknown error 37",					// ERROR37
			"File name too long",				// ENAMETOOLONG
			"No locks available",				// ENOLCK
			"Function not implemented",			// ENOSYS
			"Directory not empty",				// ENOTEMPTY
			"Illegal byte sequence",			// EILSEQ
			"Unknown error 43",					// ERROR43
			"Unknown error 44",					// ERROR44
			"Unknown error 45",					// ERROR45
			"Unknown error 46",					// ERROR46
			"Unknown error 47",					// ERROR47
			"Unknown error 48",					// ERROR48
			"Unknown error 49",					// ERROR49
			"Invalid packet size",				// EPACKSIZE
			"Not enough buffers left",			// EOUTOFBUFS
			"Illegal ioctl for device",			// EBADIOCTL
			"Bad mode for ioctl",				// EBADMODE
			"Would block",						// EWOULDBLOCK
			"Bad destination address",			// EBADDEST
			"Destination not reachable",		// EDSTNOTRCH
			"Already connected",				// EISCONN
			"Address in use",					// EADDRINUSE
			"Connection refused",				// ECONNREFUSED
			"Connection reset",					// ECONNRESET
			"Connection timed out",				// ETIMEDOUT
			"Urgent data present",				// EURG
			"No urgent data present",			// ENOURG
			"No connection",					// ENOTCONN
			"Already shutdown",					// ESHUTDOWN
			"No such connection",				// ENOCONN
			"Address family not supported",		// EAFNOSUPPORT
			"Protocol not supported by AF",		// EPROTONOSUPPORT
			"Protocol wrong type for socket",	// EPROTOTYPE
			"Operation in progress",			// EINPROGRESS
			"Address not available",			// EADDRNOTAVAIL
			"Connection already in progress",	// EALREADY
			"Message too long",					// EMSGSIZE
			"Socket operation on non-socket",	// ENOTSOCK
			"Protocol not available",			// ENOPROTOOPT
			"Unknown error 76",					// ERROR76
			"Unknown error 77",					// ERROR77
			"Unknown error 78",					// ERROR78
			"Unknown error 79",					// ERROR79
			"String was truncated",				// STRUNCATE
		};

		static readonly int _sys_nerr=_sys_errlist.Length;

		public static int errno=0;

		public static string strerror(ERRORNUMBER err)
		{
			return strerror((int)err);
		}

		public static string strerror(int err)
		{
			if(err>=0||err<_sys_nerr) return _sys_errlist[err];
			return "Unknown Error "+err;
		}
	}
}
