using System;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		// determine latitude from authalic latitude
		public static double[] pj_authset(double es)
		{
			const double P00=0.33333333333333333333; // 1/3
			const double P01=0.17222222222222222222; // 31/180
			const double P02=0.10257936507936507937; // 517/5040
			const double P10=0.06388888888888888888; // 23/360
			const double P11=0.06640211640211640212; // 251/3780
			const double P20=0.01677689594356261023; // 761/45360

			const int APA_SIZE=3;

			try
			{
				double t;
				double[] APA=new double[APA_SIZE];

				APA[0]=es*P00;
				t=es*es;
				APA[0]+=t*P01;
				APA[1]=t*P10;
				t*=es;
				APA[0]+=t*P02;
				APA[1]+=t*P11;
				APA[2]=t*P20;

				return APA;
			}
			catch
			{
				return null;
			}
		}

		public static double pj_authlat(double beta, double[] APA)
		{
			double t=beta+beta;
			return beta+APA[0]*Math.Sin(t)+APA[1]*Math.Sin(t+t)+APA[2]*Math.Sin(t+t+t);
		}
	}
}
