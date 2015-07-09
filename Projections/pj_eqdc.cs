using System;
using System.Text;
using Free.Ports.Proj4.FactorsDeriv;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_eqdc : PJ, IFactors
	{
		protected double phi1, phi2, n, rho, rho0, c;
		protected double[] en;
		protected bool ellips;

		public override string Name { get { return "eqdc"; } }
		public override string DescriptionName { get { return "Equidistant Conic"; } }
		public override string DescriptionType { get { return "Conic, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "lat_1= lat_2="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_1={0}", phi1*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lat_2={0}", phi2*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		public void_LP_PJ_FACTORS spc { get; protected set; }

		const double EPS10=1.0e-10;

		// sphere & ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			rho=c-(ellips?Proj.pj_mlfn(lp.phi, Math.Sin(lp.phi), Math.Cos(lp.phi), en):lp.phi);
			lp.lam*=n;
			xy.x=rho*Math.Sin(lp.lam);
			xy.y=rho0-rho*Math.Cos(lp.lam);

			return xy;
		}

		// sphere & ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y=rho0-xy.y;
			rho=Libc.hypot(xy.x, xy.y);
			if(rho!=0.0)
			{
				if(n<0.0)
				{
					rho=-rho;
					xy.x=-xy.x;
					xy.y=-xy.y;
				}
				lp.phi=c-rho;
				if(ellips) lp.phi=Proj.pj_inv_mlfn(ctx, lp.phi, es, en);
				lp.lam=Math.Atan2(xy.x, xy.y)/n;
			}
			else
			{
				lp.lam=0.0;
				lp.phi=n>0.0?Proj.HALFPI:-Proj.HALFPI;
			}

			return lp;
		}

		void fac(LP lp, ref FACTORS fac)
		{
			double sinphi=Math.Sin(lp.phi);
			double cosphi=Math.Cos(lp.phi);
			fac.code|=Factors.IS_ANAL_HK;
			fac.h=1.0;
			fac.k=n*(c-(ellips?Proj.pj_mlfn(lp.phi, sinphi, cosphi, en):lp.phi))/Proj.pj_msfn(sinphi, cosphi, es);
		}

		public override PJ Init()
		{
			phi1=Proj.pj_param_r(ctx, parameters, "lat_1");
			phi2=Proj.pj_param_r(ctx, parameters, "lat_2");

			if(Math.Abs(phi1+phi2)<EPS10) { Proj.pj_ctx_set_errno(ctx, -21); return null; }

			en=Proj.pj_enfn(es);
			if(en==null) return null;

			double sinphi=n=Math.Sin(phi1);
			double cosphi=Math.Cos(phi1);
			bool secant=Math.Abs(phi1-phi2)>=EPS10;
			ellips=(es>0.0);
			if(ellips)
			{
				double m1=Proj.pj_msfn(sinphi, cosphi, es);
				double ml1=Proj.pj_mlfn(phi1, sinphi, cosphi, en);
				if(secant)
				{ // secant cone
					sinphi=Math.Sin(phi2);
					cosphi=Math.Cos(phi2);
					n=(m1-Proj.pj_msfn(sinphi, cosphi, es))/(Proj.pj_mlfn(phi2, sinphi, cosphi, en)-ml1);
				}
				c=ml1+m1/n;
				rho0=c-Proj.pj_mlfn(phi0, Math.Sin(phi0), Math.Cos(phi0), en);
			}
			else
			{
				if(secant) n=(cosphi-Math.Cos(phi2))/(phi2-phi1);
				c=phi1+Math.Cos(phi1)/n;
				rho0=c-phi0;
			}

			inv=e_inverse;
			fwd=e_forward;
			spc=fac;

			return this;
		}
	}
}
