using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_imw_p : PJ
	{
		protected double P, Pp, Q, Qp, R_1, R_2, sphi_1, sphi_2, C2, phi_1, phi_2, lam_1;
		protected double[] en;
		protected int mode; // 0 if phi_1 and phi_2!=0; 1 if phi_1==0; -1 if phi_2==0

		public override string Name { get { return "imw_p"; } }
		public override string DescriptionName { get { return "International Map of the World Polyconic"; } }
		public override string DescriptionType { get { return "Mod. Polyconic, Ell"; } }
		public override string DescriptionParameters { get { return "lat_1= and lat_2= [lon_1=]"; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_1={0}", phi_1*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lon_1={0}", lam_1*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lat_2={0}", phi_2*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		const double TOL=1.0e-10;
		const double EPS=1.0e-10;

		int phi12(out double del, out double sig)
		{
			if(!Proj.pj_param_t(ctx, parameters, "lat_1")||!Proj.pj_param_t(ctx, parameters, "lat_2"))
			{
				del=sig=0;
				return -41;
			}

			phi_1=Proj.pj_param_r(ctx, parameters, "lat_1");
			phi_2=Proj.pj_param_r(ctx, parameters, "lat_2");
			del=0.5*(phi_2-phi_1);
			sig=0.5*(phi_2+phi_1);
			return (Math.Abs(del)<EPS||Math.Abs(sig)<EPS)?-42:0;
		}

		XY loc_for(LP lp, out double yc)
		{
			XY xy;

			yc=0;

			if(lp.phi==0)
			{
				xy.x=lp.lam;
				xy.y=0.0;
				return xy;
			}

			double sp=Math.Sin(lp.phi);
			double m=Proj.pj_mlfn(lp.phi, sp, Math.Cos(lp.phi), en);
			double xa=Pp+Qp*m;
			double ya=P+Q*m;
			double R=1.0/(Math.Tan(lp.phi)*Math.Sqrt(1.0-es*sp*sp));
			double C=Math.Sqrt(R*R-xa*xa);
			if(lp.phi<0.0) C=-C;
			C+=ya-R;

			double xb, yb;

			if(mode<0)
			{
				xb=lp.lam;
				yb=C2;
			}
			else
			{
				double t=lp.lam*sphi_2;
				xb=R_2*Math.Sin(t);
				yb=C2+R_2*(1.0-Math.Cos(t));
			}

			double xc;
			if(mode>0)
			{
				xc=lp.lam;
				yc=0.0;
			}
			else
			{
				double t=lp.lam*sphi_1;
				xc=R_1*Math.Sin(t);
				yc=R_1*(1.0-Math.Cos(t));
			}

			double D=(xb-xc)/(yb-yc);
			double B=xc+D*(C+R-yc);
			xy.x=D*Math.Sqrt(R*R*(1+D*D)-B*B);
			if(lp.phi>0) xy.x=-xy.x;
			xy.x=(B+xy.x)/(1.0+D*D);
			xy.y=Math.Sqrt(R*R-xy.x*xy.x);
			if(lp.phi>0) xy.y=-xy.y;
			xy.y+=C+R;

			return xy;
		}

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double yc;
			return loc_for(lp, out yc);
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=phi_2;
			lp.lam=xy.x/Math.Cos(lp.phi);

			XY t;
			do
			{
				double yc;
				t=loc_for(lp, out yc);
				lp.phi=((lp.phi-phi_1)*(xy.y-yc)/(t.y-yc))+phi_1;
				lp.lam=lp.lam*xy.x/t.x;
			} while(Math.Abs(t.x-xy.x)>TOL||Math.Abs(t.y-xy.y)>TOL);

			return lp;
		}

		void xy(double phi, out double x, out double y, ref double sp, ref double R)
		{
			sp=Math.Sin(phi);
			R=1.0/(Math.Tan(phi)*Math.Sqrt(1.0-es*sp*sp));
			double F=lam_1*sp;
			y=R*(1-Math.Cos(F));
			x=R*Math.Sin(F);
		}

		public override PJ Init()
		{
			en=Proj.pj_enfn(es);
			if(en==null) return null;

			double del, sig;
			int i=phi12(out del, out sig);
			if(i!=0) { Proj.pj_ctx_set_errno(ctx, i); return null; }

			if(phi_2<phi_1)
			{ // make sure P->phi_1 most southerly
				del=phi_1;
				phi_1=phi_2;
				phi_2=del;
			}

			if(Proj.pj_param_t(ctx, parameters, "lon_1")) lam_1=Proj.pj_param_r(ctx, parameters, "lon_1");
			else
			{ // use predefined based upon latitude
				sig=Math.Abs(sig*Proj.RAD_TO_DEG);
				if(sig<=60) sig=2.0;
				else if(sig<=76) sig=4.0;
				else sig=8.0;
				lam_1=sig*Proj.DEG_TO_RAD;
			}

			double x1, x2, T2, y1;
			mode=0;

			if(phi_1!=0) xy(phi_1, out x1, out y1, ref sphi_1, ref R_1);
			else
			{
				mode=1;
				y1=0.0;
				x1=lam_1;
			}

			if(phi_2!=0) xy(phi_2, out x2, out T2, ref sphi_2, ref R_2);
			else
			{
				mode=-1;
				T2=0.0;
				x2=lam_1;
			}

			double m1=Proj.pj_mlfn(phi_1, sphi_1, Math.Cos(phi_1), en);
			double m2=Proj.pj_mlfn(phi_2, sphi_2, Math.Cos(phi_2), en);
			double t=m2-m1;
			double s=x2-x1;
			double y2=Math.Sqrt(t*t-s*s)+y1;

			C2=y2-T2;
			t=1.0/t;
			P=(m2*y1-m1*y2)*t;
			Q=(y2-y1)*t;
			Pp=(m2*x1-m1*x2)*t;
			Qp=(x2-x1)*t;
			fwd=e_forward;
			inv=e_inverse;

			return this;
		}
	}
}
