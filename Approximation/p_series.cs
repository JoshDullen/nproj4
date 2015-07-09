#if USE_NOT_LIB_STUFF
using System.IO;

// print row coefficients of Tseries structure

namespace Free.Ports.Proj4.Approximation
{
	public static partial class Approximation
	{
		public static void p_series(Tseries T, TextWriter file, string fmt)
		{
			const int CUT=60;	// check length of line
			string format=" {0:"+fmt+"}%n";

			file.WriteLine("u: {0}", T.mu+1);
			for(int i=0; i<=T.mu; i++)
			{
				if(T.cu[i].m!=0)
				{
					string tmp=string.Format("{0} {1}", i, T.cu[i].m);
					file.Write(tmp);
					int L=tmp.Length;
					int n=0;
					for(int j=0; j<T.cu[i].m; j++)
					{
						L+=n;
						if(L>CUT)
						{
							file.Write("\n ");
							L=1;
						}
						tmp=string.Format(format, T.cu[i].c[j]);
						file.Write(tmp);
						n=tmp.Length;
					}
					file.WriteLine();
				}
			}

			file.WriteLine("v: {0}", T.mv+1);
			for(int i=0; i<=T.mv; i++)
			{
				if(T.cv[i].m!=0)
				{
					string tmp=string.Format("{0} {1}", i, T.cv[i].m);
					file.Write(tmp);
					int L=tmp.Length;
					int n=0;
					for(int j=0; j<T.cv[i].m; j++)
					{
						L+=n;
						if(L>CUT)
						{
							file.Write("\n ");
							L=1;
						}
						tmp=string.Format(format, T.cv[i].c[j]);
						file.Write(tmp);
						n=tmp.Length;
					}
					file.WriteLine();
				}
			}
		}
	}
}
#endif
