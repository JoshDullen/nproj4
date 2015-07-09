using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_cass : PJ
	{
		protected double m0;
		protected double[] en;

		public override string Name { get { return "cass"; } }
		public override string DescriptionName { get { return "Cassini"; } }
		public override string DescriptionType { get { return "Cyl, Sph&Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double C1=0.16666666666666666666;
		const double C2=0.00833333333333333333;
		const double C3=0.04166666666666666666;
		const double C4=0.33333333333333333333;
		const double C5=0.06666666666666666666;

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double n=Math.Sin(lp.phi);
			double c=Math.Cos(lp.phi);
			xy.y=Proj.pj_mlfn(lp.phi, n, c, en);

			n=1.0/Math.Sqrt(1.0-es*n*n);
			double tn=Math.Tan(lp.phi);
			double t=tn*tn;
			double a1=lp.lam*c;
			c*=es*c/(1-es);
			double a2=a1*a1;
			xy.x=n*a1*(1.0-a2*t*(C1-(8.0-t+8.0*c)*a2*C2));
			xy.y-=m0-n*tn*a2*(0.5+(5.0-t+6.0*c)*a2*C3);

			return xy;
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=Math.Asin(Math.Cos(lp.phi)*Math.Sin(lp.lam));
			xy.y=Math.Atan2(Math.Tan(lp.phi), Math.Cos(lp.lam))-phi0;

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double ph1=Proj.pj_inv_mlfn(ctx, m0+xy.y, es, en);
			double tn=Math.Tan(ph1);
			double t=tn*tn;
			double n=Math.Sin(ph1);
			double r=1.0/(1.0-es*n*n);
			n=Math.Sqrt(r);
			r*=(1.0-es)*n;
			double dd=xy.x/n;
			double d2=dd*dd;
			lp.phi=ph1-(n*tn/r)*d2*(0.5-(1.0+3.0*t)*d2*C3);
			lp.lam=dd*(1.0+t*d2*(-C4+(1.0+3.0*t)*d2*C5))/Math.Cos(ph1);

			return lp;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double dd=xy.y+phi0;
			lp.phi=Math.Asin(Math.Sin(dd)*Math.Cos(xy.x));
			lp.lam=Math.Atan2(Math.Tan(xy.x), Math.Cos(dd));

			return lp;
		}

		public override PJ Init()
		{
			if(es!=0)
			{
				en=Proj.pj_enfn(es);
				if(en==null) return null;
				m0=Proj.pj_mlfn(phi0, Math.Sin(phi0), Math.Cos(phi0), en);
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
