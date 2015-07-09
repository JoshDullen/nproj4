using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_loxim : PJ
	{
		protected double phi1, cosphi1, tanphi1;

		public override string Name { get { return "loxim"; } }
		public override string DescriptionName { get { return "Loximuthal"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return "lat_1="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_1={0}", phi1*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		const double EPS8=1.0e-8;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.y=lp.phi-phi1;
			if(Math.Abs(xy.y)<EPS8) xy.x=lp.lam*cosphi1;
			else
			{
				xy.x=Proj.FORTPI+0.5*lp.phi;
				if(Math.Abs(xy.x)<EPS8||Math.Abs(Math.Abs(xy.x)-Proj.HALFPI)<EPS8) xy.x=0.0;
				else xy.x=lp.lam*xy.y/Math.Log(Math.Tan(xy.x)/tanphi1);
			}

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y+phi1;
			if(Math.Abs(xy.y)<EPS8) lp.lam=xy.x/cosphi1;
			else
			{
				lp.lam=Proj.FORTPI+0.5*lp.phi;
				if(Math.Abs(lp.lam)<EPS8||Math.Abs(Math.Abs(lp.lam)-Proj.HALFPI)<EPS8) lp.lam=0.0;
				else lp.lam=xy.x*Math.Log(Math.Tan(lp.lam)/tanphi1)/xy.y;
			}

			return lp;
		}

		public override PJ Init()
		{
			phi1=Proj.pj_param_r(ctx, parameters, "lat_1");
			cosphi1=Math.Cos(phi1);
			if(cosphi1<EPS8) { Proj.pj_ctx_set_errno(ctx, -22); return null; }
			tanphi1=Math.Tan(Proj.FORTPI+0.5*phi1);
			inv=s_inverse;
			fwd=s_forward;
			es=0.0;

			return this;
		}
	}
}
