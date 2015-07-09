using System;

namespace Free.Ports.Proj4.LibcStuff
{
	public static partial class Libc
	{
		public const double HUGE_VAL=double.MaxValue;

		// Convert nptr to a double.
		// newptr contains the characters after the last one converted.
		public static double strtod(string nptr, out string newptr)
		{
			const double DBL_MIN=2.2250738585072014e-308;

			newptr=nptr; // if there was no number to convert.

			try
			{
				if(nptr==null)
				{
					errno=(int)ERRORNUMBER.EINVAL;
					return 0;
				}

				string s=nptr;
				if(s.Length==0)
				{
					errno=(int)ERRORNUMBER.EINVAL;
					return 0;
				}

				s=s.TrimStart(); // Skip white space.
				if(s.Length==0)
				{
					errno=(int)ERRORNUMBER.EINVAL;
					return 0;
				}

				// Check for a sign.
				short sign=1;
				if(s[0]=='-')
				{
					sign=-1;
					s=s.Substring(1);
				}
				else if(s[0]=='+') s=s.Substring(1);

				if(s.Length==0) return 0;

				double num=0;			// The number so far.
				bool got_dot=false;		// Found a decimal point.
				bool got_digit=false;	// Seen any digits.
				long exponent=0;		// The exponent of the number.

				while(s.Length>0)
				{
					if('0'<=s[0]&&s[0]<='9')
					{
						got_digit=true;

						// Make sure that multiplication by 10 will not overflow.
						if(num>double.MaxValue*0.1)
						{
							// The value of the digit doesn't matter, since we have already
							// gotten as many digits as can be represented in a 'double'.
							// This doesn't necessarily mean the result will overflow.
							// The exponent may reduce it to within range.
							//
							// We just need to record that there was another
							// digit so that we can multiply by 10 later.
							exponent++;
						}
						else num=num*10.0+(s[0]-'0');

						// Keep track of the number of digits after the decimal point.
						// If we just divided by 10 here, we would lose precision.
						if(got_dot) exponent--;
					}
					else if(!got_dot&&s[0]=='.') got_dot=true; // Record that we have found the decimal point.
					else break; // Any other character terminates the number.

					s=s.Substring(1);
				}

				if(!got_digit) return 0;

				if(s.Length>0&&(s[0]=='e'||s[0]=='E'||s[0]=='d'||s[0]=='D'))
				{
					// Get the exponent specified after the 'e' or 'E' or 'd' or 'D'.
					string end;
					long exp;

					char s0=s[0];

					s=s.Substring(1);
					exp=strtol(s, out end, 10);

					if(errno==(int)ERRORNUMBER.ERANGE)
					{
						// The exponent overflowed a 'int'. It is probably a safe
						// assumption that an exponent that cannot be represented by
						// a 'int' exceeds the limits of a 'double'.
						newptr=end;
						if(exp<0) return 0;
						return HUGE_VAL*sign;
					}

					// There was no exponent.
					if(end.Length==s.Length) end=s0+s;

					s=end;
					exponent+=exp;
				}

				newptr=s;

				if(num==0.0) return 0.0;

				// Multiply NUM by 10 to the EXPONENT power,
				// checking for overflow and underflow.
				if(exponent<0)
				{
					if(num<DBL_MIN*Math.Pow(10, -exponent))
					{
						errno=(int)ERRORNUMBER.ERANGE;
						return 0;
					}
				}
				else if(exponent>0)
				{
					if(num>double.MaxValue*Math.Pow(10, -exponent))
					{
						errno=(int)ERRORNUMBER.ERANGE;
						return HUGE_VAL*sign;
					}
				}

				num*=Math.Pow(10, exponent);
				return num*sign;
			}
			catch
			{
				errno=(int)ERRORNUMBER.EINVAL;
				return 0;
			}
		}
	}
}
