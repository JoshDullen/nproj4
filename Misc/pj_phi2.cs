using System;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// determine latitude angle phi-2
		public static double pj_phi2(projCtx ctx, double ts, double e)
		{
			const int N_ITER=15;

			double eccnth=0.5*e;
			double Phi=HALFPI-2.0*Math.Atan(ts);
			int i=N_ITER;
			double dphi=0;
			do
			{
				double con=e*Math.Sin(Phi);
				dphi=HALFPI-2.0*Math.Atan(ts*Math.Pow((1.0-con)/(1.0+con), eccnth))-Phi;
				Phi+=dphi;
			} while(Math.Abs(dphi)>TOL10&&(--i)!=0);

			if(i<=0) pj_ctx_set_errno(ctx, -18);

			return Phi;
		}
	}
}
