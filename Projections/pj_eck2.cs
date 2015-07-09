using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_eck2 : PJ
	{
		public override string Name { get { return "eck2"; } }
		public override string DescriptionName { get { return "Eckert II"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double FXC=0.46065886596178063902;
		const double FYC=1.44720250911653531871;

		const double C13=0.33333333333333333333;
		const double ONEEPS=1.0000001;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.y=Math.Sqrt(4.0-3.0*Math.Sin(Math.Abs(lp.phi)));
			xy.x=FXC*lp.lam*xy.y;
			xy.y=FYC*(2.0-xy.y);

			if(lp.phi<0.0) xy.y=-xy.y;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=2.0-Math.Abs(xy.y)/FYC;
			lp.lam=xy.x/(FXC*lp.phi);
			lp.phi=(4.0-lp.phi*lp.phi)*C13;

			if(Math.Abs(lp.phi)>=1.0)
			{
				if(Math.Abs(lp.phi)>ONEEPS) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				else lp.phi=lp.phi<0.0?-Proj.HALFPI:Proj.HALFPI;
			}
			else lp.phi=Math.Asin(lp.phi);

			if(xy.y<0) lp.phi=-lp.phi;

			return lp;
		}

		public override PJ Init()
		{
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}
