using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_wag3 : PJ
	{
		protected double C_x;

		public override string Name { get { return "wag3"; } }
		public override string DescriptionName { get { return "Wagner III"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return "lat_ts="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();

				double x=Proj.HALFPI/2;
				if(C_x<=0) x=Proj.HALFPI;
				else if(C_x>=1) x=0;
				else
				{
					double diff=0;
					do
					{
						double x3=x/3;
						diff=(Math.Cos(x)/Math.Cos(2*x3)-C_x)*(3*Math.Cos(4*x3)+3)/(Math.Sin(5*x3)+5*Math.Sin(x3));
						x=x+diff;
					} while(Math.Abs(diff)>EPS);
				}

				ret.AppendFormat(nc, " +lat_ts={0}", x*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		const double EPS=1.0e-10;
		const double TWOTHIRD_wag3=0.6666666666666666666667;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=C_x*lp.lam*Math.Cos(TWOTHIRD_wag3*lp.phi);
			xy.y=lp.phi;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y;
			lp.lam=xy.x/(C_x*Math.Cos(TWOTHIRD_wag3*lp.phi));

			return lp;
		}

		public override PJ Init()
		{
			double ts=Proj.pj_param_r(ctx, parameters, "lat_ts");
			C_x=Math.Cos(ts)/Math.Cos(2.0*ts/3.0);
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}
