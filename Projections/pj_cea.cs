using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_cea : PJ
	{
		protected double qp;
		protected double[] apa;

		public override string Name { get { return "cea"; } }
		public override string DescriptionName { get { return "Equal Area Cylindrical"; } }
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

		// spheroid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=k0*lp.lam;
			xy.y=0.5*Proj.pj_qsfn(Math.Sin(lp.phi), e, one_es)/k0;
			return xy;
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=k0*lp.lam;
			xy.y=Math.Sin(lp.phi)/k0;
			return xy;
		}

		// spheroid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=Proj.pj_authlat(Math.Asin(2.0*xy.y*k0/qp), apa);
			lp.lam=xy.x/k0;
			return lp;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y*=k0;
			double t=Math.Abs(xy.y);
			if(t-EPS10<=1.0)
			{
				if(t>=1.0) lp.phi=xy.y<0.0?-Proj.HALFPI:Proj.HALFPI;
				else lp.phi=Math.Asin(xy.y);
				lp.lam=xy.x/k0;
			}
			else Proj.pj_ctx_set_errno(ctx, -20);

			return lp;
		}

		public override PJ Init()
		{
			double t=0;

			if(Proj.pj_param_t(ctx, parameters, "lat_ts"))
			{
				t=Proj.pj_param_r(ctx, parameters, "lat_ts");
				k0=Math.Cos(t);
				if(k0<0.0) { Proj.pj_ctx_set_errno(ctx, -24); return null; }
			}

			if(es!=0)
			{
				t=Math.Sin(t);
				k0/=Math.Sqrt(1.0-es*t*t);
				e=Math.Sqrt(es);
				apa=Proj.pj_authset(es);
				if(apa==null) return null;
				qp=Proj.pj_qsfn(1.0, e, one_es);
				inv=e_inverse;
				fwd=e_forward;
			}
			else
			{
				inv=s_inverse;
				fwd=s_forward;
			}

			return this;
		}
	}
}
