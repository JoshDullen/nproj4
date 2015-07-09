#if USE_NOT_LIB_STUFF
using System;
using Free.Ports.Proj4.LibcStuff;

// Error message processing
namespace Free.Ports.Proj4.Approximation
{
	public partial class Approximation
	{
		public struct EMESS
		{
			public string File_name;	// input file name
			public string Prog_name;	// name of program
			public int File_line;		// approximate line read where error occured
		}

		// for emess procedure
		public static EMESS emess_dat=new EMESS();

		public static void emess(int code, string fmt, params object[] args)
		{
			// prefix program name, if given
			if(fmt!=null) Console.Error.Write("{0}\n<{1}>: ", Proj.pj_get_release(), emess_dat.Prog_name);

			// print file name and line, if given
			if(emess_dat.File_name!=null&&emess_dat.File_name.Length!=0)
			{
				Console.Error.Write("while processing file: {0}", emess_dat.File_name);
				if(emess_dat.File_line>0) Console.Error.WriteLine(", line {0}", emess_dat.File_line);
				else Console.Error.WriteLine();
			}
			else Console.Error.WriteLine();

			// if |code|==2, print errno code data
			if(code==2||code==-2) Console.Error.WriteLine("Sys errno: {0}: {1}", Libc.errno, Libc.strerror(Libc.errno));

			// post remainder of call data
			Console.Error.Write(fmt, args);
			
			// die if code positive
			if(code>0)
			{
				Console.Error.WriteLine("\nprogram abnormally terminated");
				throw new Exception(code.ToString());
			}
			else Console.Error.WriteLine();
		}
	}
}
#endif
