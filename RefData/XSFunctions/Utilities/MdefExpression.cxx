#include <XSUtil/Numerics/MathOperator.h>
#include <XSUtil/Utils/IosHolder.h>
#include <XSUtil/Utils/XSstream.h>
#include <XSUtil/Utils/XSutility.h>
#include <cctype>
#include <cmath>
#include <stack>
#include <utility>

// MdefExpression
#include <XSFunctions/Utilities/MdefExpression.h>
#include <XSFunctions/Utilities/funcType.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <XSFunctions/Utilities/XSModelFunction.h>

string MdefElementString[] = {"ENG", "ENGC", "NUM", "PARAM", "OPER", "UFUNC", "BFUNC", 
			      "LPAREN", "RPAREN", "COMMA", "XSMODEL", "CONXSMODEL"};

// Access to the list of models

extern ModelFunctionMap XSFunctionMap;

// Class MdefExpression::MdefExpressionError 

MdefExpression::MdefExpressionError::MdefExpressionError (const string& errMsg)
   :YellowAlert("\nMdefine Expression Error: ")
{
  *IosHolder::errHolder() << errMsg << std::endl;
}


// Class MdefExpression 
const string MdefExpression::s_allValidChars = string("_#,.+-/*^():abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789 \t\r\n");
MdefExpression::MathOpContainer MdefExpression::s_operatorsMap;
std::map<string,int> MdefExpression::s_precedenceMap;

MdefExpression::MdefExpression(const MdefExpression &right)
   : AbstractExpression(right),
     m_distinctParNames(right.m_distinctParNames),
     m_paramsToGet(right.m_paramsToGet),
     m_paramTokenIndex(right.m_paramTokenIndex),
     m_numericalConsts(right.m_numericalConsts),
     m_operators(right.m_operators),
     m_postfixElems(right.m_postfixElems),
     m_infixElems(right.m_infixElems),
     m_eLow(right.m_eLow),
     m_eHigh(right.m_eHigh),
     m_compType(right.m_compType)
{
   if (s_operatorsMap.empty())
      buildOperatorsMap();
}

MdefExpression::MdefExpression (std::pair<Real,Real> eLimits, const string& compType)
   : AbstractExpression(),
     m_distinctParNames(),
     m_paramsToGet(),
     m_paramTokenIndex(),
     m_numericalConsts(),
     m_operators(),
     m_postfixElems(),
     m_infixElems(),
     m_eLow(eLimits.first),
     m_eHigh(eLimits.second),
     m_compType(compType)
{
   if (s_operatorsMap.empty())
      buildOperatorsMap();
}


MdefExpression::~MdefExpression()
{
}


MdefExpression & MdefExpression::operator=(const MdefExpression &right)
{
   if (this != &right)
   {
      MdefExpression tmp(right);
      Swap(tmp);
   }
   return *this;
}


void MdefExpression::init (const string& exprString, bool removeWhitespace)
{
   AbstractExpression::init(exprString, removeWhitespace);
   convertToInfix();
   fixConvolutionComponents();
   convertToPostfix();
}

void MdefExpression::Swap (MdefExpression& right)
{
   AbstractExpression::Swap(right);
   std::swap(m_distinctParNames,right.m_distinctParNames);
   std::swap(m_numericalConsts,right.m_numericalConsts);
   std::swap(m_postfixElems,right.m_postfixElems);
   std::swap(m_infixElems,right.m_infixElems);
   std::swap(m_operators,right.m_operators);
   std::swap(m_paramsToGet,right.m_paramsToGet);
   std::swap(m_paramTokenIndex,right.m_paramTokenIndex);
   std::swap(m_eLow,right.m_eLow);
   std::swap(m_eHigh,right.m_eHigh);
   std::swap(m_compType,right.m_compType);
}

MdefExpression* MdefExpression::clone () const
{
   return new MdefExpression(*this);
}

const string& MdefExpression::allValidChars () const
{
   return s_allValidChars;
}

void MdefExpression::buildOperatorsMap ()
{
   using namespace Numerics;
   clearOperatorsMap();
   // s_operatorsMap will own the memory of these objects throughout
   // the program's lifetime.
   s_operatorsMap["+"] = new PlusOp();
   s_operatorsMap["-"] = new MinusOp();
   s_operatorsMap["*"] = new MultOp();
   s_operatorsMap["/"] = new DivideOp();
   s_operatorsMap["^"] = new PowOp();
   s_operatorsMap["max"] = new MaxOp();
   s_operatorsMap["min"] = new MinOp();
   s_operatorsMap["@"] = new UnaryMinusOp();
   s_operatorsMap["exp"] = new ExpOp();
   s_operatorsMap["sin"] = new SinOp();
   s_operatorsMap["sind"] = new SinDOp();
   s_operatorsMap["cos"] = new CosOp();
   s_operatorsMap["cosd"] = new CosDOp();
   s_operatorsMap["tan"] = new TanOp();
   s_operatorsMap["tand"] = new TanDOp();
   s_operatorsMap["log"] = new LogOp();
   s_operatorsMap["ln"] = new LnOp();
   s_operatorsMap["sqrt"] = new SqrtOp();
   s_operatorsMap["abs"] = new AbsOp();
   s_operatorsMap["int"] = new IntOp();
   s_operatorsMap["asin"] = new ASinOp();
   s_operatorsMap["acos"] = new ACosOp();
   s_operatorsMap["mean"] = new MeanOp();
   s_operatorsMap["dim"] = new DimOp();
   s_operatorsMap["smin"] = new SMinOp();
   s_operatorsMap["smax"] = new SMaxOp();

   s_precedenceMap["+"] = 0;
   s_precedenceMap["-"] = 0;
   s_precedenceMap["@"] = 0;
   s_precedenceMap["*"] = 1;
   s_precedenceMap["/"] = 1;
   s_precedenceMap["^"] = 2;
}

void MdefExpression::clearOperatorsMap ()
{
   MathOpContainer::iterator itOp = s_operatorsMap.begin();
   MathOpContainer::iterator itOpEnd = s_operatorsMap.end();
   while (itOp != itOpEnd)   
   {
      delete itOp->second;
      ++itOp;
   }
   s_operatorsMap.clear();
   s_precedenceMap.clear();
}

void MdefExpression::convertToInfix ()
{
   // Basically convert tokens enumerated by AbstractExpression's Token
   // types into MdefExpression's ElementTypes enumerators, which are more
   // useful for analyzing equations.
   std::vector<size_t> nonNumberTokens;
   findTheNumbers(nonNumberTokens,m_numericalConsts);
   int prevIdx = -1;
   int idx = -2;
   size_t numCount = 0;
   m_paramTokenIndex.clear();
   for (size_t i=0; i<nonNumberTokens.size(); ++i)
   {
      idx = static_cast<int>(nonNumberTokens[i]);
      MdefExpression::ElementType mathType = ENG; // init is irrelevant.
      // Assume any sequence of skipped tokens have been bundled into
      // one and only one number.
      if (idx > prevIdx+1)
      {
         m_infixElems.push_back(NUM);
         mathType = NUM;
         ++numCount;
      }
      const TokenType& curTok = tokenList()[idx];
      switch (curTok.type)
      {
         case WordExp:
            mathType = classifyWords(curTok.tokenString);
	    if ( mathType == PARAM ) {
	      m_paramTokenIndex.push_back(idx);
	    }
            break;
         case Lbrace:
            // Conditions for implied '*': '(' is preceded by
            // ')' OR Eng OR param name OR NUM.
            if (m_infixElems.size())
            {
	      MdefExpression::ElementType prev = m_infixElems[m_infixElems.size()-1];
               if (prev == RPAREN || prev == ENG || prev == PARAM 
                        || prev == NUM || prev == ENGC)
               {
                  m_infixElems.push_back(OPER);
                  m_operators.push_back("*");
               }
            }
            mathType = LPAREN;
            break;
         case Rbrace:
            // This may also include an implied '*', which is
            // inserted after exiting the switch block.
            mathType = RPAREN;
            break;
         case Plus:
            mathType = OPER;
            m_operators.push_back("+");
            break;
         case Minus:
            mathType = OPER;
            // Conditions for unary: First token OR preceded by
            // a *,/,(, or Comma.
            if (idx == 0)
               m_operators.push_back("@");
            else
            {
               Token preType = tokenList()[idx-1].type;
               if (preType == Lbrace || preType == Star || 
                        preType == Slash || preType == Comma)
                  m_operators.push_back("@");
               else
                  m_operators.push_back("-");
            }
            break;
         case Star:
            mathType = OPER;
            m_operators.push_back("*");
            break;
         case Slash:
            mathType = OPER;
            m_operators.push_back("/");
            break;
         case Exp:
            mathType = OPER;
            m_operators.push_back("^");
            break;
         case Comma:
            mathType = COMMA;
            break;
         default:
            {
               string errMsg("Unrecognized symbol during infix parsing: ");
               errMsg += curTok.tokenString;
               throw MdefExpressionError(errMsg);
            }
            break;
      }
      m_infixElems.push_back(mathType);
      // Conditions for implied '*' after ')': ')' is followed by
      // a WordExp of any kind.  
      if (mathType == RPAREN && static_cast<int>(tokenList().size()) > idx+1)
      {
         if (tokenList()[idx+1].type == WordExp)
         {
            m_infixElems.push_back(OPER);
            m_operators.push_back("*");
         }
      }
      prevIdx = idx;
   } // End non-NUM token loop

   // idx may be negative here, so we don't want to cast it into a size_t.
   if (idx+1 < static_cast<int>(tokenList().size()))
   {
      // Assume the last token must be a number.  (There should never
      // be more than one remaining number, a condition which should 
      // already have been filtered out by AbstractExpression parsing.)
      m_infixElems.push_back(NUM);
      ++numCount;
   }
   if (numCount != m_numericalConsts.size())
     throw MdefExpressionError("Last symbol is not a number");

   // now we need to do some rearrangement in the case of any CONXSMODEL instances. In general
   // a convolution model will be written before the things it acts on however we need to
   // change this so a convolution model will be written after the things it acts on.

   // Start 
   
   verifyInfix();

   std::ostringstream oss;
   oss << "Infix elements: ";
   for (size_t i=0; i<m_infixElems.size(); ++i)
     oss << MdefElementString[m_infixElems[i]] << " ";
   oss << std::endl << "Distinct parameter names: ";
   for (size_t i=0; i<m_distinctParNames.size(); ++i)
     oss << m_distinctParNames[i] << " ";
   oss << std::endl << "Parameters to get: ";
   for (size_t i=0; i<m_paramsToGet.size(); ++i)
     oss << m_paramsToGet[i] << " ";
   oss << std::endl << "Parameters token index: ";
   for (size_t i=0; i<m_paramTokenIndex.size(); ++i)
     oss << m_paramTokenIndex[i] << " ";
   oss << std::endl << "Numerical consts: ";
   for (size_t i=0; i<m_numericalConsts.size(); ++i)
     oss << m_numericalConsts[i] << " ";
   oss << std::endl << "Infix operators: ";
   for (size_t i=0; i<m_operators.size(); ++i)
     oss << m_operators[i] << " ";
   oss << std::endl;
   FunctionUtility::xsWrite(oss.str(), 40);
   
}

void MdefExpression::fixConvolutionComponents()
{

  size_t numElems(m_infixElems.size());
  
  // it is useful to have arrays of parameters to get, parameters token index,
  // numerical constants, and operator names on the Infix element array
  size_t opPos = 0;
  size_t parPos = 0;
  size_t numPos = 0;
  std::vector<string> opArray(numElems, " ");
  std::vector<int> parToGetArray(numElems, -1);
  std::vector<Real> numConstArray(numElems, -999.);

  for (size_t iElem=0; iElem<numElems; ++iElem) {
    if ( m_infixElems[iElem] == PARAM ) {
      parToGetArray[iElem] = m_paramsToGet[parPos];
      parPos++;
    }
    if ( m_infixElems[iElem] == NUM ) {
      numConstArray[iElem] = m_numericalConsts[numPos];
      numPos++;
    }
    if ( m_infixElems[iElem] == OPER || m_infixElems[iElem] == XSMODEL ||
	 m_infixElems[iElem] == CONXSMODEL ) {
      opArray[iElem] = m_operators[opPos];
      opPos++;
    }
  }

  // loop through Infix elements
  for (size_t iElem=0; iElem<numElems; ++iElem) {

    // check if this element is a convolution component
    if ( m_infixElems[iElem] == CONXSMODEL ) {
      size_t firstConvCompElem = iElem;
      // first find the end of the entire convolution component which will correspond to
      // parenCount being 0 since we know that firstConvCompElem+1 must be a LPAREN.
      int parenCount(1);
      size_t lastConvCompElem=firstConvCompElem+2;
      while ( parenCount != 0 && lastConvCompElem < numElems-1 ) {
	if ( m_infixElems[lastConvCompElem] == LPAREN ) parenCount++;
	if ( m_infixElems[lastConvCompElem] == RPAREN ) parenCount--;
	lastConvCompElem++;
      }
      if ( lastConvCompElem != numElems-1 ) lastConvCompElem--;

      // if lastConvCompElem is not a * then this convolution component has already been
      // processed so is already in the correct place

      if ( opArray[lastConvCompElem+1] == "*" ) {

	// loop through elements until we find a "+" or a "-" when the parenCount is zero
	// or the parenCount goes negative starting at lastConvCompElem+2 since we know
	// that lastConvCompElem with be "*"
	size_t firstElem=lastConvCompElem+2;
	size_t lastElem=firstElem;
	parenCount = 0;
	while ( !(parenCount == 0 && (opArray[lastElem] == "+" || opArray[lastElem] == "-")) &&
		parenCount >= 0 && lastElem < numElems-1 ) {
	  if ( m_infixElems[lastElem] == LPAREN ) parenCount++;
	  if ( m_infixElems[lastElem] == RPAREN ) parenCount--;
	  lastElem++;
	}
	if ( parenCount < 0 ) lastElem--;
	if ( lastElem != numElems-1 ) lastElem--;

	// create temporary versions of the arrays
	std::vector<MdefExpression::ElementType> tinfixElems(numElems);
	std::vector<string> topArray(numElems);
	std::vector<int> tparToGetArray(numElems);
	std::vector<Real> tnumConstArray(numElems);

	// fill temporary versions of arrays from the original arrays using:
	//    0 -> firstConvCompElem-1
	//    firstElem -> lastElem
	//    lastConvCompElem+1
	//    firstConvCompElem -> lastConvCompElem
	//    lastElem + 1 -> end

	size_t jpt(0);
	if ( firstConvCompElem > 0 ) {
	  for (size_t jElem=0; jElem<=firstConvCompElem-1; jElem++) {
	    tinfixElems[jpt] = m_infixElems[jElem];
	    topArray[jpt] = opArray[jElem];
	    tparToGetArray[jpt] = parToGetArray[jElem];
	    tnumConstArray[jpt] = numConstArray[jElem];
	    jpt++;
	  }
	}
	for (size_t jElem=firstElem; jElem<=lastElem; jElem++) {
	  tinfixElems[jpt] = m_infixElems[jElem];
	  topArray[jpt] = opArray[jElem];
	  tparToGetArray[jpt] = parToGetArray[jElem];
	  tnumConstArray[jpt] = numConstArray[jElem];
	  jpt++;
	}
	tinfixElems[jpt] = m_infixElems[lastConvCompElem+1];
	topArray[jpt] = opArray[lastConvCompElem+1];
	tparToGetArray[jpt] = parToGetArray[lastConvCompElem+1];
	tnumConstArray[jpt] = numConstArray[lastConvCompElem+1];
	jpt++;
	for (size_t jElem=firstConvCompElem; jElem<=lastConvCompElem; jElem++) {
	  tinfixElems[jpt] = m_infixElems[jElem];
	  topArray[jpt] = opArray[jElem];
	  tparToGetArray[jpt] = parToGetArray[jElem];
	  tnumConstArray[jpt] = numConstArray[jElem];
	  jpt++;
	}
	for (size_t jElem=lastElem+1; jElem<numElems; jElem++) {
	  tinfixElems[jpt] = m_infixElems[jElem];
	  topArray[jpt] = opArray[jElem];
	  tparToGetArray[jpt] = parToGetArray[jElem];
	  tnumConstArray[jpt] = numConstArray[jElem];
	  jpt++;
	}

	// now put the contents of the temporary arrays back into the original arrays
	for (size_t jElem=0; jElem<numElems; jElem++) {
	  m_infixElems[jElem] = tinfixElems[jElem];
	  opArray[jElem] = topArray[jElem];
	  parToGetArray[jElem] = tparToGetArray[jElem];
	  numConstArray[jElem] = tnumConstArray[jElem];
	}
	// end if lastElem >= firstElem
      }
      // end if CONXSMODEL  
    }
    // end iElem loop
  }

  // reconstruct the m_paramsToGet, m_numericalConsts, and m_operators arrays
  // also update m_paramTokenIndex
  opPos = 0;
  parPos = 0;
  numPos = 0;
  for (size_t iElem=0; iElem<numElems; iElem++) {
    if ( m_infixElems[iElem] == PARAM ) {
      m_paramsToGet[parPos] = parToGetArray[iElem];
      m_paramTokenIndex[parPos] = iElem;
      parPos++;
    }
    if ( m_infixElems[iElem] == NUM ) {
      m_numericalConsts[numPos] = numConstArray[iElem];
      numPos++;
    }
    if ( m_infixElems[iElem] == OPER || m_infixElems[iElem] == XSMODEL ||
	 m_infixElems[iElem] == CONXSMODEL ) {
      m_operators[opPos] = opArray[iElem];
      opPos++;
    }
  }

   std::ostringstream oss;
   oss << "fixConvolutionComponents: ";
   for (size_t i=0; i<m_infixElems.size(); ++i)
     oss << MdefElementString[m_infixElems[i]] << " ";
   oss << std::endl << "Distinct parameter names: ";
   for (size_t i=0; i<m_distinctParNames.size(); ++i)
     oss << m_distinctParNames[i] << " ";
   oss << std::endl << "Parameters to get: ";
   for (size_t i=0; i<m_paramsToGet.size(); ++i)
     oss << m_paramsToGet[i] << " ";
   oss << std::endl << "Parameters token index: ";
   for (size_t i=0; i<m_paramTokenIndex.size(); ++i)
     oss << m_paramTokenIndex[i] << " ";
   oss << std::endl << "Numerical consts: ";
   for (size_t i=0; i<m_numericalConsts.size(); ++i)
     oss << m_numericalConsts[i] << " ";
   oss << std::endl << "Operators: ";
   for (size_t i=0; i<m_operators.size(); ++i)
     oss << m_operators[i] << " ";
   oss << std::endl;
   FunctionUtility::xsWrite(oss.str(), 40);

}


void MdefExpression::convertToPostfix ()
{
   // This will fill in the m_postfixElems vector and also rearrange 
   // m_operators into the order they will be called in postfix.
   // The vectors for pars and numbers need no such reordering since 
   // they are called in the same sequence for infix and postfix.
   using namespace std;
   using Numerics::MathOperator;
   // The int in the pair below is for storing operator precedence number.
   stack<pair<int, string> > opStack;
   vector<string> tmpOperators;
   const size_t nElems = m_infixElems.size();
   size_t opPos = 0;
   for (size_t i=0; i<nElems; ++i) {
     MdefExpression::ElementType curType = m_infixElems[i];
     switch (curType) {
     case OPER:
       {
	 string curOp(m_operators[opPos]);
	 int prec = s_precedenceMap.find(curOp)->second;
	 if (!opStack.empty() && curOp != "^") {
	   pair<int,string> topOp = opStack.top();
	   int testPrec = topOp.first;
	   while (testPrec >= prec) {
	     m_postfixElems.push_back(OPER);
	     tmpOperators.push_back(topOp.second);
	     opStack.pop();
	     if (opStack.empty())
	       testPrec = -999; // cause immediate exit from loop
	     else {
	       topOp = opStack.top();
	       testPrec = topOp.first;
	     }
	   }               
	 }
	 opStack.push(make_pair(prec, curOp));
	 ++opPos;
       }
       break;
     case UFUNC:
     case BFUNC:
     case XSMODEL:
     case CONXSMODEL:
       // Unary, binary, and xspec function calls will be handled the same way.  
       // Treat them as an LPAREN that also happens to have an actual operator 
       // associated with it. This means giving it a precedence of -1 so that 
       // RPAREN thinks it's a matching LPAREN, BUT also give it the actual 
       // function name rather than the blank string.
       opStack.push(make_pair(-1,m_operators[opPos]));
       ++opPos;
       // Now skip the actual LPAREN that we know is following this.
       ++i;
       break;
     case LPAREN:
       // If in here, this is a stand-alone parenthesis not associated with a 
       // function call.  LPAREN is not an actual operator, but we still need 
       // to store it as a placeholder in the opStack.  Give it a dummy 
       // precedence that's lower than any real operator, and a blank string.
       opStack.push(make_pair(-1," "));
       break;
     case RPAREN:
       // Pop the operators stack until we come to the first left parenthesis,
       // which may or may not be tied to a function call.
       {
	 int testPrec = 0;
	 do {
	   pair<int,string> topOp = opStack.top();
	   testPrec = topOp.first;
	   if (topOp.second != " ") {
	     // This could be an LPAREN tied to a function call, if testPrec
	     // = -1.  From this point forward, BFUNC and UFUNC need only be 
	     // classified as an OPER.
	     m_postfixElems.push_back(OPER);
	     tmpOperators.push_back(topOp.second);
	   }
	   opStack.pop();  
	 } while (testPrec != -1 && !opStack.empty());
       }
       break;
     case COMMA:
       // Similar to RPAREN case, but do NOT pop the function call LPAREN.
       {
	 pair<int,string> topOp = opStack.top();
	 int testPrec = topOp.first;
	 // We can safely assume the first LPAREN reached is the 
	 // crucial function call LPAREN.  Any intervening ones would
	 // already have been popped when their RPAREN was processed.
	 while (testPrec != -1) {
	   m_postfixElems.push_back(OPER);
	   tmpOperators.push_back(topOp.second);
	   opStack.pop();
	   topOp = opStack.top();
	   testPrec = topOp.first;
	 }
       }
       break;
     default:
       m_postfixElems.push_back(curType);
       break;
     }
   }

   // Any operators remaining on the stack must now be popped to output.
   while (!opStack.empty())
   {
      m_postfixElems.push_back(OPER);
      tmpOperators.push_back(opStack.top().second);
      opStack.pop();
   }

   // And now the costly copy...
   m_operators = tmpOperators;

   std::ostringstream oss;
   oss << "Postfix elements: ";
   for (size_t i=0; i<m_postfixElems.size(); ++i)
     oss << MdefElementString[m_postfixElems[i]] << " ";
   oss << std::endl << "Postfix operators: ";
   for (size_t i=0; i<m_operators.size(); ++i)
     oss << m_operators[i] << " ";
   oss << std::endl;
   FunctionUtility::xsWrite(oss.str(), 40);

}

MdefExpression::ElementType MdefExpression::classifyWords (const string& wordStr)
{
   // Can assume wordStr is not a number, those have been taken care of
   // already.  It must either be a function, 'e', 'E', (with optional '.'
   // if convolution) or parameter name (which must start with a letter). 

   MdefExpression::ElementType type = COMMA; // Init to something other than ENG/ENGC 
   const bool isCon = (m_compType == string("con")); 
   if (wordStr.size() == 1 && (wordStr[0] == 'e' || wordStr[0] == 'E'))
      type = isCon ? ENGC : ENG;       
   else if (isCon && (wordStr == string(".e") || wordStr == string(".E")))
      type = ENG; 

   if (type != ENG && type != ENGC)
   {  
      // Function match will be case-insensitive, but other than that it
      // must match exactly.
      string lcWord = XSutility::lowerCase(wordStr);
      MathOpContainer::const_iterator itFunc = s_operatorsMap.find(lcWord);
      if (itFunc !=  s_operatorsMap.end() && lcWord == itFunc->first)
      {
         type = (itFunc->second->nArgs() == 1) ? UFUNC : BFUNC;
         m_operators.push_back(lcWord);
      } else {
	// It's not a function so check for XSPEC model using full name resolution
	bool foundXSPECmodel(false);
	if ( wordStr.length() > 2 ) {
	  try {
	    ComponentInfo nameFound = XSModelFunction::fullMatch(wordStr);
	    if ( nameFound.type() == string("con") ) {
	      type = CONXSMODEL;
	    } else {
	      type = XSMODEL;
	    }
	    m_operators.push_back(nameFound.name());
	    foundXSPECmodel = true;
	  } catch(...) {
	  }
	}
	if ( !foundXSPECmodel ) {
	  // At this point, must assume wordStr is a parameter name.
	  if (!(isalpha(wordStr[0]) || wordStr[0] == '_' || 
		wordStr.find(":") != string::npos))
	    {
	      string errMsg("Illegal parameter name: ");
	      errMsg += wordStr;
	      throw MdefExpressionError(errMsg);
	    }
	  // Determine if this name has already been entered, and if so,
	  // what is its index.   This is O(N^2) (ugh), but can't 
	  // imagine N > 100.
	  size_t idx=0;
	  std::vector<string>::const_iterator itName = m_distinctParNames.begin();
	  std::vector<string>::const_iterator itNameEnd = m_distinctParNames.end();
	  while (itName != itNameEnd)
	    {
	      if (*itName == wordStr)
		break;
	      ++itName, ++idx;
	    }
	  if (itName == itNameEnd) {
	    m_distinctParNames.push_back(wordStr);
	  }
	  m_paramsToGet.push_back(idx);
	  type = PARAM;
	}
      } 
   } // end if not ENG/ENGC
   return type;
}

void MdefExpression::verifyInfix () const
{
   // AbstractExpression has already verified general-syntax Token
   // ordering, which is really most of the work.  This will perform 
   // some additional checking that is specific to equation requirements:
   // - All function and xspec model names must immediately be followed by a '('.
   // - Commas can only exist in binary function or xspec model calls
   // - 1 and only 1 comma must exist in (the top level of) a binary
   //   function call
   // - Comma cannot appear inside any pair of parentheses internal to
   //   the pair specifying the binary function call to which it belongs.

   // This function needs something added to check that the number of commas
   // for an xspec model corresponds to the number of parameters for that model


   const size_t nElems = m_infixElems.size();
   size_t funcCounter = 0; // Needed only for output messages.
   // The first rule is easy enough to check with a linear search.
   // While doing this, count up all the commas and binary function
   // calls.  That these be equal is necessary but not sufficient
   // for rules 2 and 3.  But once rule 3 is verified independently,
   // this becomes sufficient for proving rule 2.  Rule 4 is verified
   // internally in verifyBFunc.
   size_t commaCount = 0;
   size_t binaryCount = 0;
   size_t xsModCommas = 0;
   for (size_t i=0; i<nElems; ++i)
   {
      MdefExpression::ElementType curType = m_infixElems[i];
      if (curType == UFUNC || curType == BFUNC)
      {
         if (i == nElems-1 || m_infixElems[i+1] != LPAREN)
         {
            string errMsg("A '(' must follow the call to: ");
            errMsg += m_operators[funcCounter];
            throw MdefExpressionError(errMsg);
         }
         ++funcCounter;
         if (curType == BFUNC)
            ++binaryCount;
      }
      else if (curType == XSMODEL || curType == CONXSMODEL)
      {
         if (i == nElems-1 || m_infixElems[i+1] != LPAREN)
         {
            string errMsg("A '(' must follow the call to: ");
	    errMsg += m_operators[funcCounter];
            throw MdefExpressionError(errMsg);
         }
         ++funcCounter;
	 xsModCommas += verifyXSmodelFunc(i, funcCounter-1);
      }
      else if (curType == OPER)
         ++funcCounter;
      else if (curType == COMMA)
         ++commaCount;
   }
   if (commaCount > binaryCount + xsModCommas)
      throw AbstractExpressionError("Commas can only exist in binary functions and xspec models.");

   if (binaryCount)
   {
      size_t idxElem = 0;
      while (idxElem < nElems)
      {
         // When verifyBFunc returns, idxElem will be set to the index
         // of the closing parenthesis to the function for which it was
         // originally called.
         if (m_infixElems[idxElem] == BFUNC)
            verifyBFunc(&idxElem);
         ++idxElem;
      }
   }   
}

void MdefExpression::verifyBFunc (size_t* idxElem) const
{
   // This is a recursive function that verifies that one and only
   // one comma is placed in a binary function call.  It also verifies
   // that the comma is not inside nested parentheses.
   // ASSUME idxElem initially is the index of a binary function, and
   // idxElem+1 is the index of an LPAREN.  These things should 
   // already have been verified.
   // When function returns, idxElem will be set to the closing
   // parenthesis of the binary function.
   size_t parCount = 1;
   size_t commasFound = 0;
   size_t idx = *idxElem+1;
   while (parCount)
   {
      ++idx;
      // Can also assume that expression does not end with a '('.
      // Therefore it should always be safe to start 2 elements 
      // further along.
      if (idx >= m_infixElems.size())
         throw RedAlert("Programming error in MdefExpression::verifyBFunc()");

      MdefExpression::ElementType curType = m_infixElems[idx];
      if (curType == COMMA)
      {
         ++commasFound;
         if (commasFound > 1)
            throw MdefExpressionError("Binary function called with more than 2 arguments");
         // This is the test which verifies rule 4 stated in verifyInfix.   
         if (parCount != 1)
            throw MdefExpressionError("Misplaced comma in binary function call.");
      }
      else if (curType == BFUNC)
      {
         verifyBFunc(&idx);
      }
      else if (curType == LPAREN)
      {
         ++parCount;
      }
      else if (curType == RPAREN)
      {
         --parCount;
      }
   }
   if (!commasFound)
   {
      throw MdefExpressionError("Binary function called with fewer than 2 arguments");
   }
   *idxElem = idx;
}

size_t MdefExpression::verifyXSmodelFunc(const size_t idxElem, const size_t opPos) const
{
  // Test that the xspec model function is used correctly in the expression
  // There must be one fewer commas than parameters in the model

  // Find the xspec model name and the number of parameters
  string opName = m_operators[opPos];
  size_t nParams = XSModelFunction::numberParameters(opName);

  // Loop forward from the current position in m_infixElems until we find
  // a matching RPAREN to the first LPAREN while counting the number of commas
  // at that level of parens.
  size_t commaCounter = 0;
  int parenLevel = 1;
  bool done = false;
  size_t iElem = idxElem+1;
  size_t nElem = m_infixElems.size();
  while ( !done ) {
    iElem++;
    if ( m_infixElems[iElem] == LPAREN ) parenLevel++;
    if ( m_infixElems[iElem] == RPAREN ) parenLevel--;
    if ( m_infixElems[iElem] == COMMA && parenLevel == 1 ) commaCounter++;
    if ( parenLevel == 0 || iElem == nElem-1 ) done = true;
  }
  if ( commaCounter != nParams-1 ) {
    std::ostringstream oss;
    oss << opName << " has " << commaCounter+1 << " arguments when " << nParams << " were expected ";
    throw MdefExpressionError(oss.str());
  }
  return commaCounter;
}

void MdefExpression::evaluate (const RealArray& energies, const RealArray& parameters, RealArray& flux, RealArray& fluxErr) const
{
  using Numerics::MathOperator;

  if (energies.size() < 2)
     throw MdefExpressionError("Energy array must be at least size 2");
  if (m_compType == string("con"))
     convolveEvaluate(energies, parameters, flux, fluxErr);
  else
  {
     const size_t nBins = energies.size() - 1;

     RealArray avgEngs(nBins);
     RealArray binWidths(nBins);
     for (size_t i=0; i<nBins; ++i)
     {
        avgEngs[i] = (energies[i+1]+energies[i])/2.0;
        binWidths[i] = fabs(energies[i+1]-energies[i]);
     }

     std::stack<RealArray> resultsStack;
     int numPos = 0;
     int parPos = 0;
     int opPos = 0;

     // m_postfixElems will contain the following enum values:
     //    ENG, NUM, PARAM, OPER
     for (size_t iElem=0; iElem<m_postfixElems.size(); iElem++)
     {
        MdefExpression::ElementType curType = m_postfixElems[iElem];

        // ENGC should never get in here, but if it does just treat
        // it like ENG.
        if (curType == ENG || curType == ENGC)
        {
	  resultsStack.push(avgEngs);
        }
        else if (curType == NUM)
        {
	  resultsStack.push(RealArray(m_numericalConsts[numPos],nBins));
          numPos++;
        }
        else if (curType == PARAM)
        {
           Real parVal = parameters[m_paramsToGet[parPos]];
	   resultsStack.push(RealArray(parVal,nBins));
           parPos++;
        }
        else if (curType == OPER)
        {
	   string opName = m_operators[opPos];
	   
	   MathOpContainer::const_iterator itFunc = s_operatorsMap.find(opName);
	   ModelFunctionMap::const_iterator itXSmod = XSFunctionMap.find(opName);
	   if ( itFunc != s_operatorsMap.end() && opName == itFunc->first) {
	     // case that OPER is a math operator
	     const MathOperator& mathFunc = *(itFunc->second);
	     const size_t nArgs = mathFunc.nArgs();
	     if (nArgs == 1) {
	       if (resultsStack.empty()) {
		 throw RedAlert("Trying to access empty stack in MdefExpression::evaluate()");
	       }
	       RealArray& top = resultsStack.top();
	       mathFunc(top);
	     } else if (nArgs == 2) {
	       // Note that second array is not a reference.
	       if (resultsStack.size() < 2)
		 throw RedAlert("Programmer error: Too few args in MdefExpression::evaluate() stack");
	       RealArray second = resultsStack.top();
	       resultsStack.pop();
	       RealArray& first = resultsStack.top();
	       mathFunc(first, second);
	     }
	   } else if ( itXSmod != XSFunctionMap.end() && opName == itXSmod->first) {
	     // case that OPER is an xspec model
	     const XSModelFunction& modFunc = *(itXSmod->second);
	     // get the number of parameters
	     size_t nParams = XSModelFunction::numberParameters(opName);
	     RealArray params(nParams);
	     // parameters will be popped off the stack in reverse order
	     for (size_t iparam=0; iparam<nParams; iparam++) {
	       params[nParams-iparam-1] = resultsStack.top()[0];
	       resultsStack.pop();
	     }
	     // find the type of this component
	     const ComponentInfo compInfo = XSModelFunction::fullMatch(opName);
	     const string m_compType = compInfo.type(); 
	     // if the component is convolution then need to pass in the current flux
	     // and flux error otherwise use new arrays
	     RealArray modFlux, modFluxErr;
	     if ( m_compType == string("con") ) {
	       if ( resultsStack.size() > 0 ) {
		 modFlux = resultsStack.top();
		 modFlux *= binWidths;
		 resultsStack.pop();
	       } else {
		 throw YellowAlert("Attempt to use a convolution component with nothing to operate on.");
	       }
	     }
	     int spectrumNumber(0);
	     string initString(" ");

	     modFunc(energies, params, spectrumNumber, modFlux, modFluxErr, initString);
	     // if the component type is add, con or mix then push the resulting
	     // flux divided by the bin width onto the stack else just push the
	     // result flux
	     if ( m_compType == string("con") || m_compType == string("add") ||
		  m_compType == string("mix") ) {
	       resultsStack.push(modFlux/binWidths);
	     } else {
	       resultsStack.push(modFlux);
	     }
	     // if this is a convolution model then want to jump the next operator
	     // and next postfixElem because the "*" is not necessary
	     if ( m_compType == string("con") ) {
	       opPos++;
	       iElem++;
	     }
	   }
	   opPos++;
	}
     } // end m_postfixElems loop

     if (resultsStack.size() != 1)
        throw RedAlert("Programmer error: MdefExpression::evaluate() stack should be of size 1 at end.");
     if (flux.size() != nBins)
        flux.resize(nBins);
     flux = resultsStack.top();

     if (m_compType != string("mul") && m_compType != string("pileup") )
     {
        // Integrate over bin, assume val is constant across bin.
        flux *= binWidths;
     }
  }
}

void MdefExpression::convolveEvaluate (const RealArray& energies, const RealArray& parameters, RealArray& flux, RealArray& fluxErr) const
{
   using Numerics::MathOperator;

   const size_t nBins = energies.size() - 1;
   if (flux.size() != nBins)
      throw RedAlert("Flux array size mismatch in mdef convolve function.");

   RealArray avgEngs(nBins);
   RealArray binWidths(nBins);
   for (size_t i=0; i<nBins; ++i)
   {
      avgEngs[i] = (energies[i+1]+energies[i])/2.0;
      binWidths[i] = fabs(energies[i+1]-energies[i]);
   }

   // To avoid building the numerical constant and parameter arrays
   // nBin times, do it once and store outside the loop (unlike the
   // standard evaluate function).
   const size_t nConsts = m_numericalConsts.size();
   const size_t nParsToGet = m_paramsToGet.size();
   std::vector<RealArray> constArrays(nConsts);
   std::vector<RealArray> paramArrays(nParsToGet);
   for (size_t i=0; i<nConsts; ++i)
   {
      constArrays[i].resize(nBins, m_numericalConsts[i]);
   }
   for (size_t i=0; i<nParsToGet; ++i)
   {
      paramArrays[i].resize(nBins, parameters[m_paramsToGet[i]]);
   }

   RealArray convFlux(0.0,nBins);
   for (size_t iBin=0; iBin<nBins; ++iBin)
   {
     // If input and convolved fluxes are row vectors [....], then
     // convEngs is the matrix convEngs_ji = avgEngs_i - avgEngs_j.
     RealArray convEngs(avgEngs[iBin] - avgEngs);
     std::stack<RealArray> resultsStack;
     size_t numPos = 0;
     size_t parPos = 0;
     size_t opPos = 0;

     // m_postfixElems will contain the following enum values:
     //    ENG, ENGC, NUM, PARAM, OPER
     for (size_t iElem=0; iElem<m_postfixElems.size(); ++iElem)
     {
        MdefExpression::ElementType curType = m_postfixElems[iElem];
        if (curType == ENG)
        {
           resultsStack.push(avgEngs);
        }
        else if (curType == ENGC)
        {
           resultsStack.push(convEngs);
        }
        else if (curType == NUM)
        {
           resultsStack.push(constArrays[numPos]);
           ++numPos;
        }
        else if (curType == PARAM)
        {
           resultsStack.push(paramArrays[parPos]);
           ++parPos;
        }
        else if (curType == OPER)
        {
	   string opName = m_operators[opPos];
	   MathOpContainer::const_iterator itFunc = s_operatorsMap.find(opName);
	   ModelFunctionMap::const_iterator itXSmod = XSFunctionMap.find(opName);
	   if ( itFunc != s_operatorsMap.end() && opName == itFunc->first) {
	     // case that OPER is a math operator
	     const MathOperator& mathFunc = *(itFunc->second);
	     const size_t nArgs = mathFunc.nArgs();
	     if (nArgs == 1) {
	       if (resultsStack.empty())
                 throw RedAlert("Trying to access empty stack in MdefExpression::convolveEvaluate()");
	       RealArray& top = resultsStack.top();
	       mathFunc(top);
	     } else if (nArgs == 2) {
	       // Note that second array is not a reference.
	       if (resultsStack.size() < 2)
                 throw RedAlert("Programmer error: Too few args in MdefExpression::convolveEvaluate() stack");
	       RealArray second = resultsStack.top();
	       resultsStack.pop();
	       RealArray& first = resultsStack.top();
	       mathFunc(first, second);
	     }
	   } else if ( itXSmod != XSFunctionMap.end() && opName == itXSmod->first) {
	     // case that OPER is an xspec model
	     const XSModelFunction& modFunc = *(itXSmod->second);
	     // get the number of parameters
	     //	     size_t nParams = modFunc.numberParameters();
	     // current workaround for the fact that we can't get the number of
	     // parameters is to define a function by eg xsfunc(param1,param2,param3,3)
	     size_t nParams = (size_t)round(resultsStack.top()[0]);
	     resultsStack.pop();
	     RealArray params(nParams);
	     // parameters will be popped off the stack in reverse order
	     for (size_t iparam=0; iparam<nParams; iparam++) {
	       params[nParams-iparam-1] = resultsStack.top()[0];
	       resultsStack.pop();
	     }
	     RealArray modFlux, modFluxErr;
	     if ( m_compType == string("con") ) {
	       if ( resultsStack.size() > 0 ) {
		 modFlux = resultsStack.top();
		 modFlux *= binWidths;
		 resultsStack.pop();
	       } else {
		 throw YellowAlert("Attempt to use a convolution component with nothing to operate on.");
	       }
	     }
	     int spectrumNumber(0);
	     string initString(" ");
	     modFunc(energies, params, spectrumNumber, modFlux, modFluxErr, initString);
	     // if the component type is add, con or mix then push the resulting
	     // flux divided by the bin width onto the stack else just push the
	     // result flux
	     if ( m_compType == string("con") || m_compType == string("add") ||
		  m_compType == string("mix") ) {
	       resultsStack.push(modFlux/binWidths);
	     } else {
	       resultsStack.push(modFlux);
	     }
	     // if this is a convolution model then want to jump the next operator
	     // and next postfixElem because the "*" is not necessary
	     if ( m_compType == string("con") ) {
	       ++opPos;
	       iElem++;
	     }
	   }
	   ++opPos;
        }
     } // end m_postfixElems loop

     if (resultsStack.size() != 1)
        throw RedAlert("Programmer error: MdefExpression::convolveEvaluate() stack should be of size 1 at end.");

     // fact is a column vector
     RealArray& fact = resultsStack.top();
     fact *= binWidths[iBin];
     // Now multiply row and col vectors for new flux.
     convFlux[iBin] = (flux*fact).sum();        
   } // end iBins loop

   flux = convFlux;
}

// Additional Declarations
