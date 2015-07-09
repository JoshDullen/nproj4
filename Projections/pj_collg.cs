using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_collg : PJ
	{
		public override string Name { get { return "collg"; } }
		public override string DescriptionName { get { return "Collignon"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double FXC=1.12837916709551257390;
		const double FYC=1.77245385090551602729;
		const double ONEEPS=1.0000001;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			if((xy.y=1.0-Math.Sin(lp.phi))<=0.0) xy.y=0.0;
			else xy.y=Math.Sqrt(xy.y);

			xy.x=FXC*lp.lam*xy.y;
			xy.y=FYC*(1.0-xy.y);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y/FYC-1.0;
			if(Math.Abs(lp.phi=1.0-lp.phi*lp.phi)<1.0) lp.phi=Math.Asin(lp.phi);
			else if(Math.Abs(lp.phi)>ONEEPS) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
			else lp.phi=lp.phi<0.0?-Proj.HALFPI:Proj.HALFPI;
			if((lp.lam=1.0-Math.Sin(lp.phi))<=0.0) lp.lam=0.0;
			else lp.lam=xy.x/(FXC*Math.Sqrt(lp.lam));
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
