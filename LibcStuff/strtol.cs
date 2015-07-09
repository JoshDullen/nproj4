namespace Free.Ports.Proj4.LibcStuff
{
	public static partial class Libc
	{
		// If base is 0 the base is determined by the presence of a leading
		// zero, indicating octal or a leading "0x" or "0X", indicating hexadecimal.
		// newptr contains the characters after the last one converted.
		public static int strtol(string nptr, out string newptr, int @base)
		{
			newptr=nptr; // if there was no number to convert.

			try
			{
				if(@base<0||@base==1||@base>36)
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
				bool negative=false;
				if(s[0]=='-')
				{
					negative=true;
					s=s.Substring(1);
				}
				else if(s[0]=='+') s=s.Substring(1);

				if(s.Length==0)
				{
					errno=(int)ERRORNUMBER.EINVAL;
					return 0;
				}

				// Recognize number prefix and if BASE is zero, figure it out ourselves.
				char s1='\0';
				if(s[0]=='0')
				{
					if(s.Length==1) return 0;
					if((@base==0||@base==16)&&(s[1]=='X'||s[1]=='x'))
					{
						s1=s[1];
						s=s.Substring(2);
						@base=16;
					}
					else if(@base==0) @base=8;
				}
				else if(@base==0) @base=10;

				uint cutoff=uint.MaxValue/(uint)@base;
				uint cutlim=uint.MaxValue%(uint)@base;

				uint i=0;

				// Save the pointer so we can check later if anything happened.
				string save=s;

				while(s.Length>0)
				{
					uint c=s[0];

					if(c>='0'&&c<='9') c-='0';
					else if(c>='a'&&c<='z') c=c-'a'+10;
					else if(c>='A'&&c<='Z') c=c-'A'+10;
					else break;

					if((int)c>=@base) break;

					// Check for overflow.
					if(i>cutoff||(i==cutoff&&c>cutlim))
					{
						errno=(int)ERRORNUMBER.ERANGE;
						if(negative) return int.MinValue;
						return int.MaxValue;
					}
					else
					{
						i*=(uint)@base;
						i+=c;
					}

					s=s.Substring(1);
				}

				// Check if anything actually happened.
				if(s.Length==save.Length)
				{
					// We must handle a special case here: the base is 0 or 16 and the
					// first two characters are '0' and 'x', but the rest are no
					// hexadecimal digits. This is no error case. We return 0 and
					// newptr points to the 'x'.
					if(s1=='x'||s1=='X') newptr=s1+save;

					return 0;
				}

				// Store in newptr the characters
				// past the last character we converted.
				newptr=s;

				// Check for a value that is within the range of 'uint', but outside the range of 'int'.
				// and Return the result of the appropriate sign.
				if(negative)
				{
					if(i<=int.MaxValue) return -(int)i;
					if((i-1)==int.MaxValue) return int.MinValue;
					errno=(int)ERRORNUMBER.ERANGE;
					return int.MinValue;
				}

				if(i>int.MaxValue)
				{
					errno=(int)ERRORNUMBER.ERANGE;
					return int.MaxValue;
				}
				return (int)i;
			}
			catch
			{
				errno=(int)ERRORNUMBER.EINVAL;
				return 0;
			}
		}
	}
}
