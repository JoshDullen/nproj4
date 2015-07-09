using System;
using Free.Ports.Proj4.ComplexPoly;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	abstract class PJ_mod_ster : PJ
	{
		public override string DescriptionType { get { return "Azi(mod)"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		protected COMPLEX[] zcoeff;
		protected double cchio, schio;
		protected int n;

		const double EPSLN=1.0e-10;

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double sinlon=Math.Sin(lp.lam);
			double coslon=Math.Cos(lp.lam);
			double esphi=e*Math.Sin(lp.phi);
			double chi=2.0*Math.Atan(Math.Tan((Proj.HALFPI+lp.phi)*0.5)*Math.Pow((1.0-esphi)/(1.0+esphi), e*0.5))-Proj.HALFPI;
			double schi=Math.Sin(chi);
			double cchi=Math.Cos(chi);
			double s=2.0/(1.0+schio*schi+cchio*cchi*coslon);

			COMPLEX p;
			p.r=s*cchi*sinlon;
			p.i=s*(cchio*schi-schio*cchi*coslon);
			p=p.pj_zpoly1(zcoeff, n);
			xy.x=p.r;
			xy.y=p.i;

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			int nn;

			COMPLEX p;
			p.r=xy.x;
			p.i=xy.y;
			for(nn=20; nn>0; nn--)
			{
				COMPLEX fpxy;
				COMPLEX fxy=p.pj_zpolyd1(zcoeff, n, out fpxy);

				fxy.r-=xy.x;
				fxy.i-=xy.y;
				double den=fpxy.r*fpxy.r+fpxy.i*fpxy.i;

				COMPLEX dp;
				dp.r=-(fxy.r*fpxy.r+fxy.i*fpxy.i)/den;
				dp.i=-(fxy.i*fpxy.r-fxy.r*fpxy.i)/den;

				p.r+=dp.r;
				p.i+=dp.i;
				if((Math.Abs(dp.r)+Math.Abs(dp.i))<=EPSLN) break;
			}

			double rh=0, sinz=0, cosz=0, phi=0;
			if(nn!=0)
			{
				rh=Libc.hypot(p.r, p.i);
				double z=2.0*Math.Atan(0.5*rh);
				sinz=Math.Sin(z);
				cosz=Math.Cos(z);
				lp.lam=lam0;
				if(Math.Abs(rh)<=EPSLN)
				{
					lp.phi=phi0;
					return lp;
				}

				double chi=Proj.aasin(ctx, cosz*schio+p.i*sinz*cchio/rh);
				phi=chi;
				for(nn=20; nn>0; nn--)
				{
					double esphi=e*Math.Sin(phi);
					double dphi=2.0*Math.Atan(Math.Tan((Proj.HALFPI+chi)*0.5)*Math.Pow((1.0+esphi)/(1.0-esphi), e*0.5))-Proj.HALFPI-phi;
					phi+=dphi;
					if(Math.Abs(dphi)<=EPSLN) break;
				}
			}

			if(nn!=0)
			{
				lp.phi=phi;
				lp.lam=Math.Atan2(p.r*sinz, rh*cchio*cosz-p.i*schio*sinz);
			}
			else lp.lam=lp.phi=Libc.HUGE_VAL;

			return lp;
		}

		// general initialization
		protected PJ setup()
		{
			double chio;

			if(es!=0)
			{
				double esphi=e*Math.Sin(phi0);
				chio=2.0*Math.Atan(Math.Tan((Proj.HALFPI+phi0)*0.5)*Math.Pow((1.0-esphi)/(1.0+esphi), e*0.5))-Proj.HALFPI;
			}
			else chio=phi0;

			schio=Math.Sin(chio);
			cchio=Math.Cos(chio);
			inv=e_inverse;
			fwd=e_forward;

			return this;
		}
	}

	class PJ_mil_os : PJ_mod_ster
	{
		public override string Name { get { return "mil_os"; } }
		public override string DescriptionName { get { return "Miller Oblated Stereographic"; } }

		public override PJ Init()
		{
			// Miller Oblated Stereographic
			COMPLEX[] AB=new COMPLEX[]
			{
				new COMPLEX(0.924500, 0.0),
				new COMPLEX(0.0, 0.0),
				new COMPLEX(0.019430, 0.0)
			};

			n=2;
			lam0=Proj.DEG_TO_RAD*20.0;
			phi0=Proj.DEG_TO_RAD*18.0;
			zcoeff=AB;
			es=0.0;

			return setup();
		}
	}

	class PJ_lee_os : PJ_mod_ster
	{
		public override string Name { get { return "lee_os"; } }
		public override string DescriptionName { get { return "Lee Oblated Stereographic"; } }

		public override PJ Init()
		{
			// Lee Oblated Stereographic
			COMPLEX[] AB=new COMPLEX[]
			{
				new COMPLEX(0.721316, 0.0),
				new COMPLEX(0.0, 0.0),
				new COMPLEX(-0.0088162, -0.00617325)
			};

			n=2;
			lam0=Proj.DEG_TO_RAD*-165.0;
			phi0=Proj.DEG_TO_RAD*-10.0;
			zcoeff=AB;
			es=0.0;

			return setup();
		}
	}

	class PJ_gs48 : PJ_mod_ster
	{
		public override string Name { get { return "gs48"; } }
		public override string DescriptionName { get { return "Mod. Stererographics of 48 U.S."; } }

		public override PJ Init()
		{
			// 48 United States
			COMPLEX[] AB=new COMPLEX[]
			{
				new COMPLEX(0.98879, 0.0),
				new COMPLEX(0.0, 0.0),
				new COMPLEX(-0.050909, 0.0),
				new COMPLEX(0.0, 0.0),
				new COMPLEX(0.075528, 0.0)
			};

			n=4;
			lam0=Proj.DEG_TO_RAD*-96.0;
			phi0=Proj.DEG_TO_RAD*-39.0;
			zcoeff=AB;
			es=0.0;
			a=6370997.0;

			return setup();
		}
	}

	class PJ_alsk : PJ_mod_ster
	{
		public override string Name { get { return "alsk"; } }
		public override string DescriptionName { get { return "Mod. Stererographics of Alaska"; } }

		public override PJ Init()
		{
			// Alaska ellipsoid
			COMPLEX[] ABe=new COMPLEX[]
			{
				new COMPLEX(0.9945303, 0.0),
				new COMPLEX(0.0052083, -0.0027404),
				new COMPLEX(0.0072721, 0.0048181),
				new COMPLEX(-0.0151089, -0.1932526),
				new COMPLEX(0.0642675, -0.1381226),
				new COMPLEX(0.3582802, -0.2884586)
			};

			// Alaska sphere
			COMPLEX[] ABs=new COMPLEX[]
			{
				new COMPLEX(0.9972523, 0.0),
				new COMPLEX(0.0052513, -0.0041175),
				new COMPLEX(0.0074606, 0.0048125),
				new COMPLEX(-0.0153783, -0.1968253),
				new COMPLEX(0.0636871, -0.1408027),
				new COMPLEX(0.3660976, -0.2937382)
			};

			n=5;
			lam0=Proj.DEG_TO_RAD*-152.0;
			phi0=Proj.DEG_TO_RAD*64.0;

			if(es!=0)
			{ // fixed ellipsoid/sphere
				zcoeff=ABe;
				a=6378206.4;
				es=0.00676866;
				e=Math.Sqrt(es);
			}
			else
			{
				zcoeff=ABs;
				a=6370997.0;
			}

			return setup();
		}
	}

	class PJ_gs50 : PJ_mod_ster
	{
		public override string Name { get { return "gs50"; } }
		public override string DescriptionName { get { return "Mod. Stererographics of 50 U.S."; } }

		public override PJ Init()
		{
			// GS50 ellipsoid
			COMPLEX[] ABe=new COMPLEX[]
			{
				new COMPLEX(0.9827497, 0.0),
				new COMPLEX(0.0210669, 0.0053804),
				new COMPLEX(-0.1031415, -0.0571664),
				new COMPLEX(-0.0323337, -0.0322847),
				new COMPLEX(0.0502303, 0.1211983),
				new COMPLEX(0.0251805, 0.0895678),
				new COMPLEX(-0.0012315, -0.1416121),
				new COMPLEX(0.0072202, -0.1317091),
				new COMPLEX(-0.0194029, 0.0759677),
				new COMPLEX(-0.0210072, 0.0834037)
			};

			// GS50 sphere
			COMPLEX[] ABs=new COMPLEX[]
			{
				new COMPLEX(0.9842990, 0.0),
				new COMPLEX(0.0211642, 0.0037608),
				new COMPLEX(-0.1036018, -0.0575102),
				new COMPLEX(-0.0329095, -0.0320119),
				new COMPLEX(0.0499471, 0.1223335),
				new COMPLEX(0.0260460, 0.0899805),
				new COMPLEX(0.0007388, -0.1435792),
				new COMPLEX(0.0075848, -0.1334108),
				new COMPLEX(-0.0216473, 0.0776645),
				new COMPLEX(-0.0225161, 0.0853673)
			};

			n=9;
			lam0=Proj.DEG_TO_RAD*-120.0;
			phi0=Proj.DEG_TO_RAD*45.0;

			if(es!=0)
			{ // fixed ellipsoid/sphere
				zcoeff=ABe;
				a=6378206.4;
				es=0.00676866;
				e=Math.Sqrt(es);
			}
			else
			{
				zcoeff=ABs;
				a=6370997.0;
			}

			return setup();
		}
	}
}
