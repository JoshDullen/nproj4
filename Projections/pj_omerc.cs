using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_omerc : PJ
	{
		//*
		protected double A, B, E, AB, ArB, BrA, rB, singam, cosgam, sinrot, cosrot, v_pole_n, v_pole_s, u_0;
		protected bool no_rot;

		public override string Name { get { return "omerc"; } }
		public override string DescriptionName { get { return "Oblique Mercator"; } }
		public override string DescriptionType { get { return "Cyl, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "no_rot alpha= [gamma=] [no_off] lonc= or lon_1= lat_1= lon_2= lat_2="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				if(no_rot) ret.Append(" +no_rot");

				double D=1.0, F=1.0;
				if(Math.Abs(phi0)>EPS)
				{
					double sinph0=Math.Sin(phi0);
					double cosph0=Math.Cos(phi0);
					D=Math.Sqrt(1.0+es*cosph0*cosph0*cosph0*cosph0/one_es)*Math.Sqrt(one_es)/(cosph0*Math.Sqrt(1.0-es*sinph0*sinph0));

					F=D*D-1.0;
					if(F<=0.0) F=0.0;
					else
					{
						F=Math.Sqrt(F);
						if(phi0<0.0) F=-F;
					}
				}

				double gamma0=Math.Atan2(singam, cosgam);
				ret.AppendFormat(nc, " +alpha={0}", Math.Asin(Math.Sin(gamma0)*D)*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lonc={0}", (lam0+Math.Asin(0.5*(F-1.0/F)*Math.Tan(gamma0))/B)*Proj.RAD_TO_DEG);

				if(u_0==0) ret.Append(" +no_off");

				return ret.ToString();
			}
		}

		const double TOL=1.0e-7;
		const double EPS=1.0e-10;

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double u, v;
			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)>EPS)
			{
				double Q=E/Math.Pow(Proj.pj_tsfn(lp.phi, Math.Sin(lp.phi), e), B);
				double temp=1.0/Q;
				double S=0.5*(Q-temp);
				double T=0.5*(Q+temp);
				double V=Math.Sin(B*lp.lam);
				double U=(S*singam-V*cosgam)/T;
				if(Math.Abs(Math.Abs(U)-1.0)<EPS) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
				v=0.5*ArB*Math.Log((1.0-U)/(1.0+U));
				temp=Math.Cos(B*lp.lam);
				if(Math.Abs(temp)<TOL) u=AB*lp.lam;
				else u=ArB*Math.Atan2(S*cosgam+V*singam, temp);
			}
			else
			{
				v=lp.phi>0?v_pole_n:v_pole_s;
				u=ArB*lp.phi;
			}
			if(no_rot)
			{
				xy.x=u;
				xy.y=v;
			}
			else
			{
				u-=u_0;
				xy.x=v*cosrot+u*sinrot;
				xy.y=u*cosrot-v*sinrot;
			}
			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double u, v;
			if(no_rot)
			{
				v=xy.y;
				u=xy.x;
			}
			else
			{
				v=xy.x*cosrot-xy.y*sinrot;
				u=xy.y*cosrot+xy.x*sinrot+u_0;
			}
			double Qp=Math.Exp(-BrA*v);
			double Sp=0.5*(Qp-1.0/Qp);
			double Tp=0.5*(Qp+1.0/Qp);
			double Vp=Math.Sin(BrA*u);
			double Up=(Vp*cosgam+Sp*singam)/Tp;
			if(Math.Abs(Math.Abs(Up)-1.0)<EPS)
			{
				lp.lam=0.0;
				lp.phi=Up<0.0?-Proj.HALFPI:Proj.HALFPI;
			}
			else
			{
				lp.phi=E/Math.Sqrt((1.0+Up)/(1.0-Up));
				if((lp.phi=Proj.pj_phi2(ctx, Math.Pow(lp.phi, 1.0/B), e))==Libc.HUGE_VAL) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				lp.lam=-rB*Math.Atan2((Sp*cosgam-Vp*singam), Math.Cos(BrA*u));
			}
			return lp;
		}

		public override PJ Init()
		{
			double alpha_c=0, gamma=0;
			no_rot=Proj.pj_param_t(ctx, parameters, "no_rot");
			bool alp=Proj.pj_param_t(ctx, parameters, "alpha");
			if(alp) alpha_c=Proj.pj_param_r(ctx, parameters, "alpha");
			bool gam=Proj.pj_param_t(ctx, parameters, "gamma");
			if(gam) gamma=Proj.pj_param_r(ctx, parameters, "gamma");

			bool no_off=false;
			double lamc=0, lam1=0, phi1=0, lam2=0, phi2=0;
			if(alp||gam)
			{
				lamc=Proj.pj_param_r(ctx, parameters, "lonc");

				no_off=Proj.pj_param_t(ctx, parameters, "no_off")||// For libproj4 compatability
					Proj.pj_param_t(ctx, parameters, "no_uoff"); // for backward compatibility
			}
			else
			{
				lam1=Proj.pj_param_r(ctx, parameters, "lon_1");
				phi1=Proj.pj_param_r(ctx, parameters, "lat_1");
				lam2=Proj.pj_param_r(ctx, parameters, "lon_2");
				phi2=Proj.pj_param_r(ctx, parameters, "lat_2");
				double con=Math.Abs(phi1);
				if(Math.Abs(phi1-phi2)<=TOL||con<=TOL||Math.Abs(con-Proj.HALFPI)<=TOL||
					Math.Abs(Math.Abs(phi0)-Proj.HALFPI)<=TOL||Math.Abs(Math.Abs(phi2)-Proj.HALFPI)<=TOL) { Proj.pj_ctx_set_errno(ctx, -33); return null; }
			}

			double com=Math.Sqrt(one_es);

			double D, F;
			if(Math.Abs(phi0)>EPS)
			{
				double sinph0=Math.Sin(phi0);
				double cosph0=Math.Cos(phi0);
				double con=1.0-es*sinph0*sinph0;
				B=cosph0*cosph0;
				B=Math.Sqrt(1.0+es*B*B/one_es);
				A=B*k0*com/con;
				D=B*com/(cosph0*Math.Sqrt(con));
				F=D*D-1.0;
				if(F<=0.0) F=0.0;
				else
				{
					F=Math.Sqrt(F);
					if(phi0<0.0) F=-F;
				}
				E=F+=D;
				E*=Math.Pow(Proj.pj_tsfn(phi0, sinph0, e), B);
			}
			else
			{
				B=1.0/com;
				A=k0;
				E=D=F=1.0;
			}

			double gamma0;
			if(alp||gam)
			{
				if(alp)
				{
					gamma0=Math.Asin(Math.Sin(alpha_c)/D);
					if(!gam) gamma=alpha_c;
				}
				else
				{
					gamma0=gamma;
					alpha_c=Math.Asin(D*Math.Sin(gamma0));
				}
				double con=Math.Abs(alpha_c);
				if(con<=TOL||Math.Abs(con-Proj.PI)<=TOL||Math.Abs(Math.Abs(phi0)-Proj.HALFPI)<=TOL) { Proj.pj_ctx_set_errno(ctx, -32); return null; }
				lam0=lamc-Math.Asin(0.5*(F-1.0/F)*Math.Tan(gamma0))/B;
			}
			else
			{
				double H=Math.Pow(Proj.pj_tsfn(phi1, Math.Sin(phi1), e), B);
				double L=Math.Pow(Proj.pj_tsfn(phi2, Math.Sin(phi2), e), B);
				F=E/H;
				double p=(L-H)/(L+H);
				double J=E*E;
				J=(J-L*H)/(J+L*H);
				double con=lam1-lam2;
				if(con<-Proj.PI) lam2-=Proj.TWOPI;
				else if(con>Proj.PI) lam2+=Proj.TWOPI;
				lam0=Proj.adjlon(0.5*(lam1+lam2)-Math.Atan(J*Math.Tan(0.5*B*(lam1-lam2))/p)/B);
				gamma0=Math.Atan(2.0*Math.Sin(B*Proj.adjlon(lam1-lam0))/(F-1.0/F));
				gamma=alpha_c=Math.Asin(D*Math.Sin(gamma0));
			}

			singam=Math.Sin(gamma0);
			cosgam=Math.Cos(gamma0);
			sinrot=Math.Sin(gamma);
			cosrot=Math.Cos(gamma);
			rB=1.0/B;
			ArB=A*rB;
			BrA=1.0/ArB;
			AB=A*B;
			if(no_off) u_0=0;
			else
			{
				u_0=Math.Abs(ArB*Math.Atan2(Math.Sqrt(D*D-1.0), Math.Cos(alpha_c)));
				if(phi0<0.0) u_0=-u_0;
			}
			F=0.5*gamma0;
			v_pole_n=ArB*Math.Log(Math.Tan(Proj.FORTPI-F));
			v_pole_s=ArB*Math.Log(Math.Tan(Proj.FORTPI+F));
			inv=e_inverse;
			fwd=e_forward;

			return this;
		}

		/*/

		protected double alpha, lamc, lam1, phi1, lam2, phi2, Gamma, al, bl, el, singam, cosgam, sinrot, cosrot, u_0;
		protected bool ellips, rot;

		public override string Name { get { return "omerc"; } }
		public override string DescriptionName { get { return "Oblique Mercator"; } }
		public override string DescriptionType { get { return "Cyl, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "no_rot rot_conv no_uoff and alpha= lonc= or lon_1= lat_1= lon_2= lat_2="; } }
		public override bool Invertable { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				if(!rot) ret.Append(" +no_rot");

				double D=1.0, F=1.0;
				if(Math.Abs(phi0)>EPS)
				{
					double sinph0=Math.Sin(phi0);
					double cosph0=Math.Cos(phi0);
					D=Math.Sqrt(1.0+es*cosph0*cosph0*cosph0*cosph0/one_es)*Math.Sqrt(one_es)/(cosph0*Math.Sqrt(1.0-es*sinph0*sinph0));

					F=D*D-1.0;
					if(F<=0.0) F=0.0;
					else
					{
						F=Math.Sqrt(F);
						if(phi0<0.0) F=-F;
					}
				}

				double gamma0=Math.Atan2(singam, cosgam);
				double sinAlpha=Math.Sin(gamma0)*D;
				ret.AppendFormat(nc, " +alpha={0}", Math.Asin(sinAlpha)*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lonc={0}", (lam0+Math.Asin(0.5*(F-1.0/F)*Math.Tan(gamma0))/bl)*Proj.RAD_TO_DEG);

				if(Math.Abs(sinrot-sinAlpha)>TOL) ret.Append(" +rot_conv");

				if(u_0==0) ret.Append(" +no_uoff");

				return ret.ToString();
			}
		}

		const double TOL=1.0e-7;
		const double EPS=1.0e-10;

		// ellipsoid & spheroid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double ul, us;
			double vl=Math.Sin(bl*lp.lam);
			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<=EPS)
			{
				ul=lp.phi<0.0?-singam:singam;
				us=al*lp.phi/bl;
			}
			else
			{
				double q=el/(ellips?Math.Pow(Proj.pj_tsfn(lp.phi, Math.Sin(lp.phi), e), bl):Math.Tan(0.5*(Proj.HALFPI-lp.phi)));
				double s=0.5*(q-1.0/q);
				ul=2.0*(s*singam-vl*cosgam)/(q+1.0/q);
				double con=Math.Cos(bl*lp.lam);
				if(Math.Abs(con)>=TOL)
				{
					us=al*Math.Atan((s*cosgam+vl*singam)/con)/bl;
					if(con<0.0) us+=Proj.PI*al/bl;
				}
				else us=al*bl*lp.lam;
			}

			if(Math.Abs(Math.Abs(ul)-1.0)<=EPS) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			double vs=0*al*Math.Log((1.0-ul)/(1.0+ul))/bl;
			us-=u_0;
			if(!rot)
			{
				xy.x=us;
				xy.y=vs;
			}
			else
			{
				xy.x=vs*cosrot+us*sinrot;
				xy.y=us*cosrot-vs*sinrot;
			}

			return xy;
		}

		// ellipsoid & spheroid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double us, vs;

			if(!rot)
			{
				us=xy.x;
				vs=xy.y;
			}
			else
			{
				vs=xy.x*cosrot-xy.y*sinrot;
				us=xy.y*cosrot+xy.x*sinrot;
			}

			us+=u_0;
			double q=Math.Exp(-bl*vs/al);

			double s=0.5*(q-1.0/q);
			double vl=Math.Sin(bl*us/al);
			double ul=2.0*(vl*cosgam+s*singam)/(q+1.0/q);
			if(Math.Abs(Math.Abs(ul)-1.0)<EPS)
			{
				lp.lam=0.0;
				lp.phi=ul<0.0?-Proj.HALFPI:Proj.HALFPI;
			}
			else
			{
				lp.phi=el/Math.Sqrt((1.0+ul)/(1.0-ul));
				if(ellips)
				{
					lp.phi=Proj.pj_phi2(ctx, Math.Pow(lp.phi, 1.0/bl), e);
					if(lp.phi==Libc.HUGE_VAL) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
				}
				else lp.phi=Proj.HALFPI-2.0*Math.Atan(lp.phi);

				lp.lam=-Math.Atan2((s*cosgam-vl*singam), Math.Cos(bl*us/al))/bl;
			}

			return lp;
		}

		public override PJ Init()
		{
			double con=0;

			rot=Proj.pj_param_b(ctx, parameters, "no_rot")==false;
			bool azi=Proj.pj_param_t(ctx, parameters, "alpha");
			if(azi)
			{
				lamc=Proj.pj_param_r(ctx, parameters, "lonc");
				alpha=Proj.pj_param_r(ctx, parameters, "alpha");
				if(Math.Abs(alpha)<=TOL||Math.Abs(Math.Abs(phi0)-Proj.HALFPI)<=TOL||Math.Abs(Math.Abs(alpha)-Proj.HALFPI)<=TOL)
				{ Proj.pj_ctx_set_errno(ctx, -32); return null; }
			}
			else
			{
				lam1=Proj.pj_param_r(ctx, parameters, "lon_1");
				phi1=Proj.pj_param_r(ctx, parameters, "lat_1");
				lam2=Proj.pj_param_r(ctx, parameters, "lon_2");
				phi2=Proj.pj_param_r(ctx, parameters, "lat_2");
				con=Math.Abs(phi1);
				if(Math.Abs(phi1-phi2)<=TOL||con<=TOL||Math.Abs(con-Proj.HALFPI)<=TOL||
					Math.Abs(Math.Abs(phi0)-Proj.HALFPI)<=TOL||Math.Abs(Math.Abs(phi2)-Proj.HALFPI)<=TOL)
				{ Proj.pj_ctx_set_errno(ctx, -33); return null; }
			}

			double d, f;

			ellips=es>0.0;
			double com=ellips?Math.Sqrt(one_es):1.0;
			if(Math.Abs(phi0)>EPS)
			{
				double sinph0=Math.Sin(phi0);
				double cosph0=Math.Cos(phi0);
				if(ellips)
				{
					con=1.0-es*sinph0*sinph0;
					bl=cosph0*cosph0;
					bl=Math.Sqrt(1.0+es*bl*bl/one_es);
					al=bl*k0*com/con;
					d=bl*com/(cosph0*Math.Sqrt(con));
				}
				else
				{
					bl=1.0;
					al=k0;
					d=1.0/cosph0;
				}

				f=d*d-1.0;

				if(f<=0.0) f=0.0;
				else
				{
					f=Math.Sqrt(f);
					if(phi0<0.0) f=-f;
				}

				f+=d;
				el=f;
				if(ellips) el*=Math.Pow(Proj.pj_tsfn(phi0, sinph0, e), bl);
				else el*=Math.Tan(0.5*(Proj.HALFPI-phi0));
			}
			else
			{
				bl=1.0/com;
				al=k0;
				el=d=f=1.0;
			}

			if(azi)
			{
				Gamma=Math.Asin(Math.Sin(alpha)/d);
				lam0=lamc-Math.Asin((0.5*(f-1.0/f))*Math.Tan(Gamma))/bl;
			}
			else
			{
				double h, l;

				if(ellips)
				{
					h=Math.Pow(Proj.pj_tsfn(phi1, Math.Sin(phi1), e), bl);
					l=Math.Pow(Proj.pj_tsfn(phi2, Math.Sin(phi2), e), bl);
				}
				else
				{
					h=Math.Tan(0.5*(Proj.HALFPI-phi1));
					l=Math.Tan(0.5*(Proj.HALFPI-phi2));
				}

				f=el/h;
				double p=(l-h)/(l+h);
				double j=el*el;
				j=(j-l*h)/(j+l*h);

				con=lam1-lam2;
				if(con<-Proj.PI) lam2-=Proj.TWOPI;
				else if(con>Proj.PI) lam2+=Proj.TWOPI;

				lam0=Proj.adjlon(0.5*(lam1+lam2)-Math.Atan(j*Math.Tan(0.5*bl*(lam1-lam2))/p)/bl);
				Gamma=Math.Atan(2.0*Math.Sin(bl*Proj.adjlon(lam1-lam0))/(f-1.0/f));
				alpha=Math.Asin(d*Math.Sin(Gamma));
			}

			singam=Math.Sin(Gamma);
			cosgam=Math.Cos(Gamma);
			f=Proj.pj_param_b(ctx, parameters, "rot_conv")?Gamma:alpha;
			sinrot=Math.Sin(f);
			cosrot=Math.Cos(f);
			u_0=Proj.pj_param_b(ctx, parameters, "no_uoff")?0.0:Math.Abs(al*Math.Atan(Math.Sqrt(d*d-1.0)/cosrot)/bl);
			if(phi0<0.0) u_0=-u_0;
			inv=e_inverse;
			fwd=e_forward;

			return this;
		}

		/**/
	}
}
