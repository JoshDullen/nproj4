using System;

// meridinal distance for ellipsoid and inverse
//	8th degree - accurate to < 1e-5 meters when used in conjuction
//		with typical major axis values.
//	Inverse determines phi to EPS (1e-11) radians, about 1e-6 seconds.

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		public static double[] pj_enfn(double es)
		{
			const double C00=1.0;
			const double C02=0.25;
			const double C04=0.046875;
			const double C06=0.01953125;
			const double C08=0.01068115234375;
			const double C22=0.75;
			const double C44=0.46875;
			const double C46=0.01302083333333333333;
			const double C48=0.00712076822916666666;
			const double C66=0.36458333333333333333;
			const double C68=0.00569661458333333333;
			const double C88=0.3076171875;

			try
			{
				double[] en=new double[5];

				double t=es*es;
				en[0]=C00-es*(C02+es*(C04+es*(C06+es*C08)));
				en[1]=es*(C22-es*(C04+es*(C06+es*C08)));
				en[2]=t*(C44-es*(C46+es*C48));
				t*=es;
				en[3]=t*(C66-es*C68);
				en[4]=t*es*C88;

				return en;
			}
			catch
			{
				return null;
			}
		}

		public static double pj_mlfn(double phi, double sphi, double cphi, double[] en)
		{
			cphi*=sphi;
			sphi*=sphi;
			return en[0]*phi-cphi*(en[1]+sphi*(en[2]+sphi*(en[3]+sphi*en[4])));
		}

		public static double pj_inv_mlfn(projCtx ctx, double arg, double es, double[] en)
		{
			const int MAX_ITER=10;

			double k=1.0/(1.0-es);
			double phi=arg;
			for(int i=MAX_ITER; i>0; --i) // rarely goes over 2 iterations
			{
				double s=Math.Sin(phi);
				double t=1.0-es*s*s;
				t=(pj_mlfn(phi, s, Math.Cos(phi), en)-arg)*t*Math.Sqrt(t)*k;
				phi-=t;
				if(Math.Abs(t)<EPS11) return phi;
			}

			pj_ctx_set_errno(ctx, -17);
			return phi;
		}
	}
}