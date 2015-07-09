using System.Collections.Generic;
using System.IO;
using System.Text;

namespace Free.Ports.Proj4.WKT
{
	public static class WKTParser
	{
		public static List<WKTNode> FromFile(string filename)
		{
			using(StreamReader sr=new StreamReader(filename)) return ParseNodes(sr);
		}

		public static List<WKTNode> FromString(string wktString)
		{
			using(StringReader sr=new StringReader(wktString)) return ParseNodes(sr);
		}

		internal static List<WKTNode> ParseNodes(TextReader tr)
		{
			bool allowComma=false;

			List<WKTNode> ret=new List<WKTNode>();

			for(; ; )
			{
				int i=tr.Peek();
				if(i<0)
				{
					if(ret.Count==0) throw new EndOfStreamException("Unexpected end of data/file, while parsing WKT.");
					return ret;
				}

				char c=(char)i;

				if(char.IsWhiteSpace(c))
				{
					tr.Read();
					continue;
				}

				if(char.IsLetter(c))
				{
					ret.Add(ParseNode(tr));
					allowComma=true;
					continue;
				}

				if(allowComma&&c==',')
				{
					tr.Read();
					allowComma=false;
					continue;
				}

				throw new InvalidDataException("Unexpected character in data/file, while parsing WKT.");
			}
		}
		
		internal static WKTNode ParseNode(TextReader tr)
		{
			StringBuilder current=new StringBuilder();
			WKTNode ret=new WKTNode();

			// Name
			bool wasWS=false;
			for(; ; )
			{
				int i=tr.Read();
				if(i<0) throw new EndOfStreamException("Unexpected end of data/file, while parsing WKT.");

				char c=(char)i;

				if(c=='[')
				{
					ret.Key=current.ToString();
					break;
				}

				if(char.IsWhiteSpace(c))
				{
					wasWS=true;
					continue;
				}

				if(wasWS||c=='"'||c==','||c==']') throw new EndOfStreamException("Unexpected end of data/file, while parsing WKT.");

				current.Append(c);
			}

			current.Clear();

			// Values
			bool allowComma=false;
			for(; ; )
			{
				int i=tr.Peek();
				if(i<0) throw new EndOfStreamException("Unexcepted end of data/file, while parsing WKT.");
				char c=(char)i;

				if(char.IsWhiteSpace(c))
				{
					tr.Read();
					continue;
				}

				if(char.IsDigit(c)||c=='+'||c=='-'||c=='.')
				{
					ret.Values.Add(ParseNumber(tr));
					allowComma=true;
					continue;
				}

				if(c=='"')
				{
					ret.Values.Add(ParseString(tr));
					allowComma=true;
					continue;
				}

				if(char.IsLetter(c))
				{
					ret.Children.Add(ParseNode(tr));
					allowComma=true;
					continue;
				}

				if(c==']')
				{
					tr.Read();
					return ret;
				}

				if(allowComma&&c==',')
				{
					allowComma=false;
					tr.Read();
					continue;
				}

				throw new InvalidDataException("Unexcepted character in data/file, while parsing WKT.");
			}
		}

		internal static string ParseNumber(TextReader tr)
		{
			StringBuilder ret=new StringBuilder();

			for(; ; )
			{
				int i=tr.Peek();
				if(i<0) throw new EndOfStreamException("Unexcepted end of data/file, while parsing WKT.");
				char c=(char)i;

				if(c==']') break;

				tr.Read();
				if(char.IsWhiteSpace(c)||c==',') break;

				if(char.IsDigit(c)||c=='+'||c=='-'||c=='.'||c=='x'||(c>='a'&&c<='f')||(c>='A'&&c<='F'))
				{
					ret.Append(c);
					continue;
				}

				throw new InvalidDataException("Unexcepted character in data/file, while parsing WKT.");
			}

			return ret.ToString();
		}

		internal static string ParseString(TextReader tr)
		{
			StringBuilder ret=new StringBuilder();

			int i=tr.Read();
			if(i<0) throw new EndOfStreamException("Unexcepted end of data/file, while parsing WKT.");
			char c=(char)i;

			if(c!='"') ret.Append(c);

			for(; ; )
			{
				i=tr.Read();
				if(i<0) throw new EndOfStreamException("Unexcepted end of data/file, while parsing WKT.");
				c=(char)i;

				if(c=='"') break;

				ret.Append(c);
			}

			return ret.ToString();
		}
	}
}
