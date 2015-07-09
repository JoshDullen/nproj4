using System;

namespace Free.Ports.Proj4.LibcStuff
{
	public class Complex
	{
		public static readonly Complex NaN=new Complex(double.NaN, double.NaN);

		public double re, im;

		public Complex(double r, double i)
		{
			re=r;
			im=i;
		}

		public static implicit operator Complex(double d)
		{
			return new Complex(d, 0);
		}

		public static Complex operator+(Complex z1, Complex z2)
		{
			return new Complex(z1.re+z2.re, z1.im+z2.im);
		}

		public static Complex operator-(Complex z1, Complex z2)
		{
			return new Complex(z1.re-z2.re, z1.im-z2.im);
		}

		public static Complex operator*(Complex z1, Complex z2)
		{
			return new Complex(z1.re*z2.re-z1.im*z2.im, z1.re*z2.im+z1.im*z2.re);
		}

		public static Complex operator/(Complex z1, Complex z2)
		{
			double c2d2=z2.re*z2.re+z2.im*z2.im;
			return new Complex((z1.re*z2.re+z1.im*z2.im)/c2d2, (z1.im*z2.re-z1.re*z2.im)/c2d2);
		}

		//public static Complex Sinh(Complex z)
		//{
		//    if(double.IsNaN(z.re)||double.IsNaN(z.im)||double.IsInfinity(z.im)) return NaN;

		//    if(z.im==0) return new Complex(Math.Sinh(z.re), 0);

		//    double coshx=Math.Cosh(z.re);
		//    if(double.IsInfinity(coshx)) return new Complex(z.re*Math.Cos(z.im), z.re*Math.Sin(z.im));
		//    return new Complex(Math.Sinh(z.re)*Math.Cos(z.im), coshx*Math.Sin(z.im));
		//}

		//public static Complex TimesI(Complex z)
		//{
		//    if(double.IsNaN(z.re)||double.IsNaN(z.im)) return NaN;
		//    if(z.re==0&&z.im==0) return z;
		//    return new Complex(-z.im, z.re);
		//}

		//public static Complex DivideI(Complex z)
		//{
		//    if(double.IsNaN(z.re)||double.IsNaN(z.im)) return NaN;
		//    if(z.re==0&&z.im==0) return z;
		//    return new Complex(z.im, -z.re);
		//}

		// The sine of this Complex.
		public static Complex Sin(Complex z)
		{
			// sin(x+iy)=sin(x)*cosh(y)+i*cos(x)*sinh(y)
			return new Complex(Math.Sin(z.re)*Math.Cosh(z.im), Math.Cos(z.re)*Math.Sinh(z.im));

			// or
			//// sin(z)=-i*sinh(i*z)
			//if(double.IsNaN(z.re)||double.IsNaN(z.im)) return NaN;
			//return DivideI(Sinh(TimesI(z)));
		}
	}
}
