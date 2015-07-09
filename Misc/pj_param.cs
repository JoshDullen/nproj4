using System.Collections.Generic;
using System.Globalization;
using Free.Ports.Proj4.LibcStuff;

// put parameters in linked list and retrieve

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		static CultureInfo nc=new CultureInfo("");

		// create parameter list entry
		internal static string pj_mkparam(string str)
		{
			try
			{
				if(str[0]=='+') return str.Substring(1);
				return str;
			}
			catch
			{
				Libc.errno=(int)ERRORNUMBER.ENOMEM;
				return null;
			}
		}

		//**********************************************************************
		//								pj_param()
		//
		//		Test for presence or get parameter value. The last
		//		character in function name is a parameter type which can take
		//		the values:
		//
		//			t - test for presence, return true/false as bool
		//			i - integer value returned as int
		//			d - simple valued real input returned as double
		//			r - degrees (DMS translation applied), returned as
		//				radians as double
		//			s - string returned as string
		//			b - test for t/T/f/F, return as bool
		//
		//**********************************************************************

		// test for presence
		internal static bool pj_param_t(projCtx ctx, List<string> pl, string name)
		{
			foreach(string p in pl)
			{
				if(p==name) return true;
				if(p.StartsWith(name+"=")) return true;
			}

			return false;
		}

		// get parameter value
		internal static bool pj_param_b(projCtx ctx, List<string> pl, string name)
		{
			if(ctx==null) ctx=pj_get_default_ctx();

			foreach(string p in pl)
			{
				if(p==name) return true;

				if(p.StartsWith(name+"="))
				{
					if(p.Length==name.Length+1) return true;
					char val=p[name.Length+1];

					switch(val)
					{
						case 'f':
						case 'F':
							return false;
						case 't':
						case 'T':
							return true;
						default:
							pj_ctx_set_errno(ctx, -8);
							return false;
					}
				}
			}

			return false;
		}

		// get parameter value
		internal static int pj_param_i(projCtx ctx, List<string> pl, string name)
		{
			if(ctx==null) ctx=pj_get_default_ctx();

			foreach(string p in pl)
			{
				if(p==name) return 0;
				if(p.StartsWith(name+"="))
				{
					if(p.Length==name.Length+1) return 0;
					string val=p.Substring(name.Length+1);

					try
					{
						return int.Parse(val);
					}
					catch
					{
						Libc.errno=(int)ERRORNUMBER.EINVAL;
						return 0;
					}
				}
			}

			return 0;
		}

		// get parameter value
		internal static double pj_param_d(projCtx ctx, List<string> pl, string name)
		{
			if(ctx==null) ctx=pj_get_default_ctx();

			foreach(string p in pl)
			{
				if(p==name) return 0;
				if(p.StartsWith(name+"="))
				{
					if(p.Length==name.Length+1) return 0;
					string val=p.Substring(name.Length+1);

					try
					{
						return double.Parse(val, nc);
					}
					catch
					{
						Libc.errno=(int)ERRORNUMBER.EINVAL;
						return 0;
					}
				}
			}

			return 0;
		}

		// get parameter value
		internal static double pj_param_r(projCtx ctx, List<string> pl, string name)
		{
			if(ctx==null) ctx=pj_get_default_ctx();

			foreach(string p in pl)
			{
				if(p==name) return 0;
				if(p.StartsWith(name+"="))
				{
					if(p.Length==name.Length+1) return 0;
					string val=p.Substring(name.Length+1);
					return dmstor_ctx(ctx, val, 0);
				}
			}

			return 0;
		}

		// get parameter value
		internal static string pj_param_s(projCtx ctx, List<string> pl, string name)
		{
			foreach(string p in pl)
			{
				if(p==name) return "";
				if(p.StartsWith(name+"="))
				{
					if(p.Length==name.Length+1) return "";
					return p.Substring(name.Length+1);
				}
			}

			return "";
		}
	}
}
