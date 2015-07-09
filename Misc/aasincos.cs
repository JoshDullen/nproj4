using System;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// arc sin, cosine, tan2 and sqrt that will NOT fail
		const double aasincos_ONE_TOL=1.00000000000001;
		
		public static double aasin(projCtx ctx, double v)
		{
			double av=Math.Abs(v);

			if(av>=1.0)
			{
				if(av>aasincos_ONE_TOL) pj_ctx_set_errno(ctx, -19);
				return v<0.0?-HALFPI:HALFPI;
			}

			return Math.Asin(v);
		}

		public static double aacos(projCtx ctx, double v)
		{
			double av=Math.Abs(v);

			if(av>=1.0)
			{
				if(av>aasincos_ONE_TOL) pj_ctx_set_errno(ctx, -19);
				return v<0.0?PI:0.0;
			}

			return Math.Acos(v);
		}

		public static double asqrt(double v)
		{
			return (v<=0)?0.0:Math.Sqrt(v);
		}

		public static double aatan2(double n, double d)
		{
			const double ATOL=1e-50;
			return (Math.Abs(n)<ATOL&&Math.Abs(d)<ATOL)?0.0:Math.Atan2(n, d);
		}
	}
}
