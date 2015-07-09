using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_mbtfpp : PJ
	{
		public override string Name { get { return "mbtfpp"; } }
		public override string DescriptionName { get { return "McBride-Thomas Flat-Polar Parabolic"; } }
		public override string DescriptionType { get { return "Cyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double CS=0.95257934441568037152;
		const double FXC=0.92582009977255146156;
		const double FYC=3.40168025708304504493;
		const double C23=0.66666666666666666666;
		const double C13=0.33333333333333333333;
		const double ONEEPS=1.0000001;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			lp.phi=Math.Asin(CS*Math.Sin(lp.phi));
			xy.x=FXC*lp.lam*(2.0*Math.Cos(C23*lp.phi)-1.0);
			xy.y=FYC*Math.Sin(C13*lp.phi);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y/FYC;
			if(Math.Abs(lp.phi)>=1.0)
			{
				if(Math.Abs(lp.phi)>ONEEPS) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				else lp.phi=(lp.phi<0.0)?-Proj.HALFPI:Proj.HALFPI;
			}
			else lp.phi=Math.Asin(lp.phi);

			lp.phi*=3.0;
			lp.lam=xy.x/(FXC*(2.0*Math.Cos(C23*lp.phi)-1.0));
			lp.phi=Math.Sin(lp.phi)/CS;
			if(Math.Abs(lp.phi)>=1.0)
			{
				if(Math.Abs(lp.phi)>ONEEPS) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				else lp.phi=(lp.phi<0.0)?-Proj.HALFPI:Proj.HALFPI;
			}
			else lp.phi=Math.Asin(lp.phi);

			return lp;
		}

		public override PJ Init()
		{
			es=0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}
