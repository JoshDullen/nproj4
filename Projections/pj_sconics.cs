using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	abstract class PJ_sconics : PJ
	{
		protected double n, rho_c, rho_0, sig, c1, c2;
		protected sconics_mode type;

		protected enum sconics_mode
		{
			EULER=0,
			MURD1=1,
			MURD2=2,
			MURD3=3,
			PCONIC=4,
			TISSOT=5,
			VITK1=6
		}
	
		public override string DescriptionType { get { return "Conic, Sph"; } }
		public override string DescriptionParameters { get { return "lat_1= and lat_2="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();

				double del=0;
				switch(type)
				{
					case sconics_mode.TISSOT: del=Math.Acos((Math.Sqrt(rho_c*rho_c-4)-rho_c)*n/2); break;
					case sconics_mode.MURD1:
						{
							double y=(rho_c-sig)*Math.Tan(sig);
							if(y>=1) throw new InvalidOperationException("Project settings invalid, or incomplete.");

							if(y<0) del=Math.PI;
							else
							{
								del=Math.Acos(y)*2;
								double diff=0;
								do
								{
									double sinx=Math.Sin(del);
									double cosx=Math.Cos(del);
									diff=(del*del*y-del*sinx)/(sinx-del*cosx);
									del=del-diff;
								} while(Math.Abs(diff)>EPS);
							}
						}
						break;
					case sconics_mode.MURD2: del=Math.Acos(rho_c*rho_c*Math.Tan(sig)*Math.Tan(sig)); break;
					case sconics_mode.MURD3:
						{
							double y=Math.Tan(sig)*(rho_c-sig);
							if(y>=1) throw new InvalidOperationException("Project settings invalid, or incomplete.");
							if(y>-0.3) del=Math.Acos(y);
							else del=2;
							double diff=0;
							do
							{
								double sinx=Math.Sin(del);
								double cosx=Math.Cos(del);
								diff=-(sinx*sinx*y-del*sinx*cosx)/(sinx*cosx-del);
								if(Math.Abs(del-diff)>=Math.PI) diff=(del-Math.PI)/2; // out of convergence zone
								del=del-diff;
							} while(Math.Abs(diff)>EPS);
						}
						break;
					case sconics_mode.EULER:
						{
							double y=n/Math.Sin(sig);
							if(y>=1||y<0) throw new InvalidOperationException("Project settings invalid, or incomplete.");

							del=Math.Acos(y)*2;
							double diff=0;
							do
							{
								double sinx=Math.Sin(del);
								double cosx=Math.Cos(del);
								diff=(del*del*y-del*sinx)/(sinx-del*cosx);
								del=del-diff;
							} while(Math.Abs(diff)>EPS);
						}
						break;
					case sconics_mode.PCONIC: del=Math.Acos(c2); break;
					case sconics_mode.VITK1:
						{
							double y=n/Math.Sin(sig);
							if(y<=1) throw new InvalidOperationException("Project settings invalid, or incomplete.");

							del=2*Math.Atan(Math.Sqrt(y-1));
							if(Math.Abs(del)>=Proj.HALFPI) del=1.5707;
							double diff=0;
							do
							{
								double sinx=Math.Sin(del);
								double cosx=Math.Cos(del);
								diff=(del*del*cosx*cosx*y-del*sinx*cosx)/(sinx*cosx-del);
								if(Math.Abs(del-diff)>=Proj.HALFPI) diff=(del-Proj.HALFPI)/2; // out of convergence zone
								del=del-diff;
							} while(Math.Abs(diff)>EPS);
						}
						break;
				}

				ret.AppendFormat(nc, " +lat_1={0}", (sig-del)*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lat_2={0}", (sig+del)*Proj.RAD_TO_DEG);

				return ret.ToString();
			}
		}

		const double EPS10=1.0e-10;
		const double EPS=1.0e-10;

		// get common factors for simple conics
		int phi12(ref double del)
		{
			if(!Proj.pj_param_t(ctx, parameters, "lat_1")||!Proj.pj_param_t(ctx, parameters, "lat_2")) return -41;

			double p1=Proj.pj_param_r(ctx, parameters, "lat_1");
			double p2=Proj.pj_param_r(ctx, parameters, "lat_2");
			del=0.5*(p2-p1);
			sig=0.5*(p2+p1);
			return (Math.Abs(del)<EPS||Math.Abs(sig)<EPS)?-42:0;
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double rho;
			switch(type)
			{
				case sconics_mode.MURD2: rho=rho_c+Math.Tan(sig-lp.phi); break;
				case sconics_mode.PCONIC: rho=c2*(c1-Math.Tan(lp.phi-sig)); break;
				default: rho=rho_c-lp.phi; break;
			}

			lp.lam*=n;
			xy.x=rho*Math.Sin(lp.lam);
			xy.y=rho_0-rho*Math.Cos(lp.lam);

			return xy;
		}

		// ellipsoid &spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y=rho_0-xy.y;
			double rho=Libc.hypot(xy.x, xy.y);
			if(n<0.0)
			{
				rho=-rho;
				xy.x=-xy.x;
				xy.y=-xy.y;
			}

			lp.lam=Math.Atan2(xy.x, xy.y)/n;

			switch(type)
			{
				case sconics_mode.PCONIC: lp.phi=Math.Atan(c1-rho/c2)+sig; break;
				case sconics_mode.MURD2: lp.phi=sig-Math.Atan(rho-rho_c); break;
				default: lp.phi=rho_c-rho; break;
			}

			return lp;
		}

		protected PJ_sconics setup()
		{
			double del=0, cs;
			int i=phi12(ref del);

			if(i!=0) { Proj.pj_ctx_set_errno(ctx, i); return null; }

			switch(type)
			{
				case sconics_mode.TISSOT:
					n=Math.Sin(sig);
					cs=Math.Cos(del);
					rho_c=n/cs+cs/n;
					rho_0=Math.Sqrt((rho_c-2*Math.Sin(phi0))/n);
					break;
				case sconics_mode.MURD1:
					rho_c=Math.Sin(del)/(del*Math.Tan(sig))+sig;
					rho_0=rho_c-phi0;
					n=Math.Sin(sig);
					break;
				case sconics_mode.MURD2:
					cs=Math.Sqrt(Math.Cos(del));
					rho_c=cs/Math.Tan(sig);
					rho_0=rho_c+Math.Tan(sig-phi0);
					n=Math.Sin(sig)*cs;
					break;
				case sconics_mode.MURD3:
					rho_c=del/(Math.Tan(sig)*Math.Tan(del))+sig;
					rho_0=rho_c-phi0;
					n=Math.Sin(sig)*Math.Sin(del)*Math.Tan(del)/(del*del);
					break;
				case sconics_mode.EULER:
					n=Math.Sin(sig)*Math.Sin(del)/del;
					del*=0.5;
					rho_c=del/(Math.Tan(del)*Math.Tan(sig))+sig;
					rho_0=rho_c-phi0;
					break;
				case sconics_mode.PCONIC:
					n=Math.Sin(sig);
					c2=Math.Cos(del);
					c1=1.0/Math.Tan(sig);
					del=phi0-sig;
					if(Math.Abs(del)-EPS10>=Proj.HALFPI) { Proj.pj_ctx_set_errno(ctx, -43); return null; }
					rho_0=c2*(c1-Math.Tan(del));
					break;
				case sconics_mode.VITK1:
					cs=Math.Tan(del);
					n=cs*Math.Sin(sig)/del;
					rho_c=del/(cs*Math.Tan(sig))+sig;
					rho_0=rho_c-phi0;
					break;
			}

			inv=s_inverse;
			fwd=s_forward;
			es=0;

			return this;
		}
	}

	class PJ_euler : PJ_sconics
	{
		public override string Name { get { return "euler"; } }
		public override string DescriptionName { get { return "Euler"; } }

		public override PJ Init()
		{
			type=sconics_mode.EULER;

			return setup();
		}
	}

	class PJ_tissot : PJ_sconics
	{
		public override string Name { get { return "tissot"; } }
		public override string DescriptionName { get { return "Tissot"; } }

		public override PJ Init()
		{
			type=sconics_mode.TISSOT;

			return setup();
		}
	}

	class PJ_murd1 : PJ_sconics
	{
		public override string Name { get { return "murd1"; } }
		public override string DescriptionName { get { return "Murdoch I"; } }

		public override PJ Init()
		{
			type=sconics_mode.MURD1;

			return setup();
		}
	}

	class PJ_murd2 : PJ_sconics
	{
		public override string Name { get { return "murd2"; } }
		public override string DescriptionName { get { return "Murdoch II"; } }

		public override PJ Init()
		{
			type=sconics_mode.MURD2;

			return setup();
		}
	}

	class PJ_murd3 : PJ_sconics
	{
		public override string Name { get { return "murd3"; } }
		public override string DescriptionName { get { return "Murdoch III"; } }

		public override PJ Init()
		{
			type=sconics_mode.MURD3;

			return setup();
		}
	}

	class PJ_pconic : PJ_sconics
	{
		public override string Name { get { return "pconic"; } }
		public override string DescriptionName { get { return "Perspective Conic"; } }

		public override PJ Init()
		{
			type=sconics_mode.PCONIC;

			return setup();
		}
	}

	class PJ_vitk1 : PJ_sconics
	{
		public override string Name { get { return "vitk1"; } }
		public override string DescriptionName { get { return "Vitkovsky I"; } }

		public override PJ Init()
		{
			type=sconics_mode.VITK1;

			return setup();
		}
	}
}
