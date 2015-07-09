using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_labrd : PJ
	{
		protected double kRg, p0s, A, C, Ca, Cb, Cc, Cd;
		protected bool rot;

		public override string Name { get { return "labrd"; } }
		public override string DescriptionName { get { return "Laborde"; } }
		public override string DescriptionType { get { return "Cyl, Sph"; } }
		public override string DescriptionParameters { get { return "azi= no_rot"; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				if(rot) ret.Append(" +rot");
				double tmp=Cb*12.0*kRg*kRg;
				if(tmp>1) tmp=1;
				else if(tmp<-1) tmp=-1;
				ret.AppendFormat(nc, " +azi={0}", Proj.RAD_TO_DEG*Math.Asin(tmp)/2);
				return ret.ToString();
			}
		}

		const double EPS=1.0e-10;

		// ellipsoid & spheroid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double V1=A*Math.Log(Math.Tan(Proj.FORTPI+0.5*lp.phi));
			double t=e*Math.Sin(lp.phi);
			double V2=0.5*e*A*Math.Log((1.0+t)/(1.0-t));
			double ps=2.0*(Math.Atan(Math.Exp(V1-V2+C))-Proj.FORTPI);
			double I1=ps-p0s;
			double cosps=Math.Cos(ps);
			double cosps2=cosps*cosps;
			double sinps=Math.Sin(ps);
			double sinps2=sinps*sinps;
			double I4=A*cosps;
			double I2=0.5*A*I4*sinps;
			double I3=I2*A*A*(5.0*cosps2-sinps2)/12.0;
			double I6=I4*A*A;
			double I5=I6*(cosps2-sinps2)/6.0;
			I6*=A*A*(5.0*cosps2*cosps2+sinps2*(sinps2-18.0*cosps2))/120.0;
			t=lp.lam*lp.lam;
			xy.x=kRg*lp.lam*(I4+t*(I5+t*I6));
			xy.y=kRg*(I1+t*(I2+t*I3));
			double x2=xy.x*xy.x;
			double y2=xy.y*xy.y;
			V1=3.0*xy.x*y2-xy.x*x2;
			V2=xy.y*y2-3.0*x2*xy.y;
			xy.x+=Ca*V1+Cb*V2;
			xy.y+=Ca*V2-Cb*V1;

			return xy;
		}

		// ellipsoid & spheroid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double x2=xy.x*xy.x;
			double y2=xy.y*xy.y;
			double V1=3.0*xy.x*y2-xy.x*x2;
			double V2=xy.y*y2-3.0*x2*xy.y;
			double V3=xy.x*(5.0*y2*y2+x2*(-10.0*y2+x2));
			double V4=xy.y*(5.0*x2*x2+y2*(-10.0*x2+y2));
			xy.x+=-Ca*V1-Cb*V2+Cc*V3+Cd*V4;
			xy.y+=Cb*V1-Ca*V2-Cd*V3+Cc*V4;
			double ps=p0s+xy.y/kRg;
			double pe=ps+phi0-p0s;

			double t;
			for(int i=20; i>0; i--)
			{
				V1=A*Math.Log(Math.Tan(Proj.FORTPI+0.5*pe));
				double tpe=e*Math.Sin(pe);
				V2=0.5*e*A*Math.Log((1.0+tpe)/(1.0-tpe));
				t=ps-2.0*(Math.Atan(Math.Exp(V1-V2+C))-Proj.FORTPI);
				pe+=t;
				if(Math.Abs(t)<EPS) break;
			}

			t=e*Math.Sin(pe);
			t=1.0-t*t;
			double Re=one_es/(t*Math.Sqrt(t));
			t=Math.Tan(ps);
			double t2=t*t;
			double s=kRg*kRg;
			double d=Re*k0*kRg;
			double I7=t/(2.0*d);
			double I8=t*(5.0+3.0*t2)/(24.0*d*s);
			d=Math.Cos(ps)*kRg*A;
			double I9=1.0/d;
			d*=s;
			double I10=(1.0+2.0*t2)/(6.0*d);
			double I11=(5.0+t2*(28.0+24.0*t2))/(120.0*d*s);
			x2=xy.x*xy.x;
			lp.phi=pe+x2*(-I7+I8*x2);
			lp.lam=xy.x*(I9+x2*(-I10+x2*I11));

			return lp;
		}

		public override PJ Init()
		{
			rot=Proj.pj_param_b(ctx, parameters, "no_rot")==false;
			double Az=Proj.pj_param_r(ctx, parameters, "azi");
			double sinp=Math.Sin(phi0);
			double t=1.0-es*sinp*sinp;
			double N=1.0/Math.Sqrt(t);
			double R=one_es*N/t;
			kRg=k0*Math.Sqrt(N*R);
			p0s=Math.Atan(Math.Sqrt(R/N)*Math.Tan(phi0));
			A=sinp/Math.Sin(p0s);
			t=e*sinp;
			C=0.5*e*A*Math.Log((1.0+t)/(1.0-t))-A*Math.Log(Math.Tan(Proj.FORTPI+0.5*phi0))+Math.Log(Math.Tan(Proj.FORTPI+0.5*p0s));
			t=Az+Az;
			Cb=1.0/(12.0*kRg*kRg);
			Ca=(1.0-Math.Cos(t))*Cb;
			Cb*=Math.Sin(t);
			Cc=3.0*(Ca*Ca-Cb*Cb);
			Cd=6.0*Ca*Cb;
			inv=e_inverse;
			fwd=e_forward;

			return this;
		}
	}
}
