using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_merc : PJ
	{
		public override string Name { get { return "merc"; } }
		public override string DescriptionName { get { return "Mercator"; } }
		public override string DescriptionType { get { return "Cyl, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "lat_ts="; } }
		public override bool Invertible { get { return true; } }

		// lat_ts is transformed to k_0 and thus don't needed anymore
		//protected override string Proj4ParameterString
		//{
		//    get
		//    {
		//        throw new DoNotImplementException();
		//    }
		//}

		const double EPS10=1.0e-10;

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<=EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
			xy.x=k0*lp.lam;
			xy.y=-k0*Math.Log(Proj.pj_tsfn(lp.phi, Math.Sin(lp.phi), e));

			return xy;
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<=EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
			xy.x=k0*lp.lam;
			xy.y=k0*Math.Log(Math.Tan(Proj.FORTPI+0.5*lp.phi));

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;
			lp.phi=Proj.pj_phi2(ctx, Math.Exp(-xy.y/k0), e);
			if(lp.phi==Libc.HUGE_VAL) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
			lp.lam=xy.x/k0;

			return lp;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=Proj.HALFPI-2.0*Math.Atan(Math.Exp(-xy.y/k0));
			lp.lam=xy.x/k0;
			return lp;
		}

		public override PJ Init()
		{
			double phits=0.0;
			bool is_phits=Proj.pj_param_t(ctx, parameters, "lat_ts");

			if(is_phits)
			{
				phits=Math.Abs(Proj.pj_param_r(ctx, parameters, "lat_ts"));
				if(phits>=Proj.HALFPI) { Proj.pj_ctx_set_errno(ctx, -24); return null; }
			}

			if(es!=0)
			{ // ellipsoid
				if(is_phits) k0=Proj.pj_msfn(Math.Sin(phits), Math.Cos(phits), es);
				inv=e_inverse;
				fwd=e_forward;
			}
			else
			{ // sphere
				if(is_phits) k0=Math.Cos(phits);
				inv=s_inverse;
				fwd=s_forward;
			}

			return this;
		}
	}
}
