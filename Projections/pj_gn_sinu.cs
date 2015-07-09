using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_gn_sinu : PJ
	{
		protected double[] en;
		protected double m, n, C_x, C_y;

		public override string Name { get { return "gn_sinu"; } }
		public override string DescriptionName { get { return "General Sinusoidal Series"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return "m= n="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +m={0}", m);
				ret.AppendFormat(nc, " +n={0}", n);
				return ret.ToString();
			}
		}

		const double EPS10=1.0e-10;
		const int MAX_ITER=8;
		const double LOOP_TOL=1.0e-7;

		// Ellipsoidal Sinusoidal only
		// ellipsoid
		protected XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double s=Math.Sin(lp.phi);
			double c=Math.Cos(lp.phi);
			xy.y=Proj.pj_mlfn(lp.phi, s, c, en);
			xy.x=lp.lam*c/Math.Sqrt(1.0-es*s*s);

			return xy;
		}

		// ellipsoid
		protected LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=Proj.pj_inv_mlfn(ctx, xy.y, es, en);
			double s=Math.Abs(lp.phi);
			if(s<Proj.HALFPI)
			{
				s=Math.Sin(lp.phi);
				lp.lam=xy.x*Math.Sqrt(1.0-es*s*s)/Math.Cos(lp.phi);
			}
			else if((s-EPS10)<Proj.HALFPI) lp.lam=0.0;
			else { Proj.pj_ctx_set_errno(ctx, -20); return lp; }

			return lp;
		}

		// General spherical sinusoidals
		// sphere
		protected XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			if(m==0) lp.phi=n!=1.0?Proj.aasin(ctx, n*Math.Sin(lp.phi)):lp.phi;
			else
			{
				double k, V;

				k=n*Math.Sin(lp.phi);
				int i=MAX_ITER;
				for(; i>0; i--)
				{
					V=(m*lp.phi+Math.Sin(lp.phi)-k)/(m+Math.Cos(lp.phi));
					lp.phi-=V;
					if(Math.Abs(V)<LOOP_TOL) break;
				}
				if(i==0) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
			}

			xy.x=C_x*lp.lam*(m+Math.Cos(lp.phi));
			xy.y=C_y*lp.phi;

			return xy;
		}

		// sphere
		protected LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y/=C_y;
			lp.phi=m!=0?Proj.aasin(ctx, (m*xy.y+Math.Sin(xy.y))/n):(n!=1.0?Proj.aasin(ctx, Math.Sin(xy.y)/n):xy.y);
			lp.lam=xy.x/(C_x*(m+Math.Cos(xy.y)));

			return lp;
		}

		// for spheres, only
		protected void setup()
		{
			es=0;
			C_y=Math.Sqrt((m+1.0)/n);
			C_x=C_y/(m+1.0);
			inv=s_inverse;
			fwd=s_forward;
		}

		public override PJ Init()
		{
			if(Proj.pj_param_t(ctx, parameters, "n")&&Proj.pj_param_t(ctx, parameters, "m"))
			{
				n=Proj.pj_param_d(ctx, parameters, "n");
				m=Proj.pj_param_d(ctx, parameters, "m");
			}
			else { Proj.pj_ctx_set_errno(ctx, -99); return null; }

			setup();

			return this;
		}
	}

	class PJ_sinu : PJ_gn_sinu
	{
		public override string Name { get { return "sinu"; } }
		public override string DescriptionName { get { return "Sinusoidal (Sanson-Flamsteed)"; } }
		public override string DescriptionType { get { return "PCyl, Sph&Ell"; } }
		public override string DescriptionParameters { get { return ""; } }

		// override PJ_gn_sinu.Proj4ParameterString as m= and n= are not needed here
		protected override string Proj4ParameterString { get { return ""; } }

		public override PJ Init()
		{
			en=Proj.pj_enfn(es);
			if(en==null) return null;

			if(es!=0)
			{
				inv=e_inverse;
				fwd=e_forward;
			}
			else
			{
				n=1.0;
				m=0.0;
				setup();
			}

			return this;
		}
	}

	class PJ_eck6 : PJ_gn_sinu
	{
		public override string Name { get { return "eck6"; } }
		public override string DescriptionName { get { return "Eckert VI"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }

		// override PJ_gn_sinu.Proj4ParameterString as m= and n= are not needed here
		protected override string Proj4ParameterString { get { return ""; } }

		public override PJ Init()
		{
			m=1.0;
			n=2.570796326794896619231321691;
			setup();

			return this;
		}
	}

	class PJ_mbtfps : PJ_gn_sinu
	{
		public override string Name { get { return "mbtfps"; } }
		public override string DescriptionName { get { return "McBryde-Thomas Flat-Polar Sinusoidal"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }

		// override PJ_gn_sinu.Proj4ParameterString as m= and n= are not needed here
		protected override string Proj4ParameterString { get { return ""; } }

		public override PJ Init()
		{
			m=0.5;
			n=1.785398163397448309615660845;
			setup();

			return this;
		}
	}
}
