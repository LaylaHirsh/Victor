// -*- C++ -*-------x-----------------------------------------------------------
//
// Description:     Abstract base class for all sorts of alignments.
//
// -----------------x-----------------------------------------------------------

#include <AlignmentBase.h>
#include <IoTools.h>

namespace Biopool
{

// CONSTRUCTORS:

AlignmentBase::AlignmentBase() : targetName(), target(), seqTemplateName(),
	seqTemplate(), startAaTarget(0)
{ }


AlignmentBase::AlignmentBase(const AlignmentBase &orig)
{
	copy(orig);
}


AlignmentBase::~AlignmentBase()
{ }


// OPERATORS:

AlignmentBase&
AlignmentBase::operator = (const AlignmentBase &orig)
{
	if (&orig != this)
		copy(orig);
	POSTCOND((orig == *this), exception);
	return *this;
}


// PREDICATES:

unsigned int
AlignmentBase::getSequenceLength(const string &seq)
{
	return (getPureSequence(seq)).length();
}


double
AlignmentBase::calculatePairwiseIdentity(const string &seq1, const string &seq2)
{
	unsigned int maxL = (seq1.length() <= seq2.length()) ? seq1.length() : seq2.length();

	if (seq1.length() != seq2.length())
		cout << "Warning: sequence lengths do not match:\n"
		     << "seq1 = " << seq1.length() << "\n"
		     << "seq2 = " << seq2.length() << "\n";

	double tmp = 0;
	int countN = 0;

	for (unsigned int i = 0; i < maxL; i++)
		if ((seq1[i] != '-') && (seq1[i] == seq2[i]))
			tmp++;

	for (unsigned int i = 0; i < maxL; i++)
		if (seq2[i] == '-')
			countN++;

//	unsigned int maxS = getSequenceLength(seq1);
//	if (maxS < getSequenceLength(seq2))
//		maxS = getSequenceLength(seq2);

//	return tmp / (getSequenceLength(seq1) - countN);
	return tmp / (seq1.length() - countN);
}


double
AlignmentBase::calculateIdentity()
{
	double tmp = 0;
	int countN = 0;

	for (unsigned int i = 0; i < target.length(); i++)
	{
		char tmpC = target[i];
		if (tmpC == '-')
			continue;

		bool add = true;
		for (unsigned int j = 0; j < seqTemplate.size(); j++)
			if (seqTemplate[j][i] == '-')
			{
				countN++;
				add = false;
				break;
			}
			else
				if (seqTemplate[j][i] != tmpC)
				{
				add = false;
				break;
				}

		if (add)
			tmp++;
	}

	return tmp / (getSequenceLength(target) - countN);
}


bool
AlignmentBase::isConserved(unsigned int p, unsigned int index) const
{
	if ((p > target.length()) || (index > this->size()))
		ERROR("Argument out of scope.", exception);

	if (index == 9999)
	{
		for (unsigned int i = 0; i < seqTemplate.size(); i++)
			if (target[p] != seqTemplate[i][p])
				return false;
	}
	else
		if (target[p] != seqTemplate[index][p])
			return false;

	return true;
}


bool
AlignmentBase::isInsertion(unsigned int p, unsigned int index) const
{
	if ((p > target.length()) || (index >= seqTemplate.size()))
		ERROR("Argument out of scope.", exception);

	if (seqTemplate[index][p] != '-')
		return false;

	return true;
}


bool
AlignmentBase::isDeletion(unsigned int p) const
{
	if (p > target.length())
		ERROR("Argument out of scope.", exception);

	if (target[p] != '-')
		return false;

	return true;
}


bool
AlignmentBase::isGap(unsigned int p, unsigned int index) const
{
	if ((p > target.length()) || (index >= seqTemplate.size()))
		ERROR("Argument out of scope.", exception);

	if ((target[p] == '-') && (seqTemplate[index][p] == '-'))
		return true;

	return false;
}


vector< vector<int> >
AlignmentBase::getMatchSubset()
{
	string targetSequence = target;
	string templateSequence = seqTemplate[0];
	string stringGap = "A-"; // gap corresponds to position 1
	int countMatchTarget = -1; // I start from minus one because the first time it enters on the loop will give position zero meaning the start of the match
	int countMatchTemplate = -1;
	vector< vector<int> > resultVecs(2);

	for (unsigned int i = 0; i < templateSequence.size(); i++)
		if (targetSequence[i] == stringGap[1])
		{
			int gap = -1;
			resultVecs[0].push_back(gap);
		}
		else
		{
			countMatchTarget++;
			resultVecs[0].push_back(countMatchTarget);
		}

	for (unsigned int i = 0; i < seqTemplate[0].size(); i++)
		if (templateSequence[i] == stringGap[1])
		{
			int gap = -1;
			resultVecs[1].push_back(gap);
		}
		else
		{
			countMatchTemplate++;
			resultVecs[1].push_back(countMatchTemplate);
		}

	return resultVecs;
}


vector<int>
AlignmentBase::shiftMatchSubset(vector<int> inputVector, int newStartPos)
{
	unsigned int len = inputVector.size();
	vector<int> newInputVector;
	int gapCount = 0;

	for (unsigned int i = 0; i < len; i++) // browse the whole vector
		if (inputVector[i] == -1)
		{
			int gap = -1;
			gapCount++;
			newInputVector.push_back(gap);
		}
		else
		{
			int valueCurrent = inputVector[i];
			int newValue = valueCurrent+ newStartPos;
			newInputVector.push_back(newValue);
		}

	return newInputVector;
}


double
AlignmentBase::matchPositionVector(vector<int> CeTarget, vector<int> CeTemplate,
	vector<int> seqTarget, vector<int> seqTemplate)
{
	double overlap = 0.00;
	int counting = 0;
	int matchCount = 0;

	// Do a cycle for knowing where the "-2" are; i think i should use a hash
	// for knowing the position and the content of that position:
	// I should take in account - what comes before (and how many times)
	// what comes after (and how many times).

	// meta code
	// preliminary loop

	// Good alignment: -2 should be always aligned with -1 for correct.

//	vector<int> placeHolderTarget;
//	for (unsigned int a = 0; a < CeTarget.size(); a++)
//		if (CeTarget == -2)
//			placeHolderTarget.push_back(a);

//	for (unsigned int b = 0; b < placeHolderTarget.size(); b += couDead)
//	{
//		int couMom = 0;
//		int couDead = 1;
//		int moment = -1;
//		int currPos = -3;
//		currPos = placeHolderTarget[b];
//		int deadZone = 0;

//		while (deadZone == 0)
//		{
//			nextPos = placeHolderTarget[b+couDead];
//			if (nextPos != (currPos + 1))
//			{
//				deadZone = 1;
//				continue;
//			}
//			currPos = nextPos;
//			couDead++;
//		}


		// At this position I know the length of the current island, namely from b
		// to b+couDead; it also means that b-1 and b+couDead+1 are the "anchors".
		// But beware, these positions can also be the end of the sequences
		// (either N terminus or C terminus).

		//       H O
		// H3-N+-C-C-N-C ...
		//       R   H

		// Check b-1 the value: we would like to know if we are on N termimus.

//		int posN = CeTarget[placeHolderTarget[b]-1]; // check which position
		// !!! VERY IMPORTANT: IF WE GO BEYOND N TERMINUS IT RETURNS SEGMENTATION FAULT! HAVE TO THINK HOW TO DEAL WITH IT !!!
//		int posC = CeTarget[placeHolderTarget[b]+couDead+1];

//		for (unsigned int c = b; c < (b + couDead); c++) // I put a very big number; don't expect any protein to reach 1000 aa (as far as I know longest sequence ever seen was around 800 aa)
//		{
//			// First thing to control is the length of the "-2" islands inside the vector
//			cou++; // increase cou;
//			moment *= -1;
//			if (CeTarget[placeHolderTarget[b]+(cou*moment)] != -2)
//			{ }
//			moment *= -1;
//		}
//	}

//	vector<int> placeHolderTemplate;
//	for (unsigned int a = 0; a < CeTemplate.size(); a++)
//		if (CeTemplate == -2)
//			placeHolderTemplate.push_back(a);

	// Once you know that position I need to know extrema to search
//	for (i = extrema1; i < extrema2; i++)
//	{
//		search for matching;
//	}

	// example 1:
	//
	// AL 1  2  3 -1 -1  4  5
	// AL 9 10 11 12 13 14 15
	//
	// CE 1  2  3 -2 -2  4  5
	// CE 9 10 11 12 13 14 15

	// example 2:
	//
	// AL 1  2 -1 -1 -1  3  4
	// AL 9 10 11 12 13 14 15
	//
	// CE 1  2 -1 -2 -2  3  4
	// CE 9 10 11 12 13 14 15

	// example 3:
	//
	// AL 1  2  3 -1 -1 -1  4
	// AL 9 10 11 12 13 14 15
	//
	// CE 1  2  3 -2 -2 -1  4
	// CE 9 10 11 12 13 14 15

	// example 4:
	//
	// AL 1  2 -1 -1 -1 -1  3
	// AL 9 10 11 12 13 14 15
	//
	// CE 1  2 -1 -2 -2 -1  3
	// CE 9 10 11 12 13 14 15


	for (unsigned int i = 0; i < seqTarget.size(); i++)
//	for (unsigned int i = 0; i < CeTarget.size(); i++)
	{
		int tmp;
		tmp = seqTarget[i];
//		tmp = CeTarget[i];

		if ((tmp == -1) || (seqTemplate[i] == -1))
			continue;
//		if ((tmp == -1) || (CeTemplate[i] == -1))
//			continue;

		matchCount++;
		for (unsigned int j = 0; j < CeTarget.size(); j++)
//		for (unsigned int j = 0; j < seqTarget.size(); j++)
			if (CeTarget[j] == tmp)
//			if (seqTarget[j] == tmp)
			{
				if (CeTemplate[j] == seqTemplate[i])
//				if (seqTemplate[j] == CeTemplate[i])
					counting++;
				break;
			}
	}

	overlap = static_cast<double>(counting) / matchCount;
	return overlap;
}


void
AlignmentBase::saveFasta(ostream &output) const
{
	saveFasta(target, targetName, output);
	for (unsigned int j = 0; j < seqTemplate.size(); j++)
		saveFasta(seqTemplate[j], seqTemplateName[j], output);
}


void
AlignmentBase::saveClustal(ostream &output) const
{
	output << "CLUSTAL\n\n";

	for (unsigned int from = 0; from < target.length(); from += 60)
	{
		saveClustal(target, this->targetName, output, from);
		for (unsigned int j = 0; j < seqTemplate.size(); j++)
			saveClustal(seqTemplate[j], this->seqTemplateName[j], output, from);
		output << "\n";
	}
}


// MODIFIERS:

void
AlignmentBase::copy(const AlignmentBase &orig)
{
	targetName = orig.targetName;
	target = orig.target;
	startAaTarget = orig.startAaTarget;

	seqTemplateName.clear();
	seqTemplate.clear();
	startAaTemplates.clear();
	for (unsigned int i = 0; i < orig.seqTemplate.size(); i++)
	{
		seqTemplateName.push_back(orig.seqTemplateName[i]);
		seqTemplate.push_back(orig.seqTemplate[i]);
		startAaTemplates.push_back(orig.startAaTemplates[i]);
	}
}


AlignmentBase*
AlignmentBase::newCopy()
{
	AlignmentBase *tmp = new AlignmentBase(*this);
	return tmp;
}


void
AlignmentBase::setTemplate(string t, string tName)
{
	if (t.length() != target.length())
		ERROR("AlignmentBase::setTemplate() Template length does not match target.", exception);
	seqTemplateName.push_back(tName);
	seqTemplate.push_back(t);
	startAaTemplates.push_back(0);
}


void
AlignmentBase::swapTemplate(unsigned int index1, unsigned int index2)
{
	if ((index1 >= this->size()) || (index2 >= this->size()))
		ERROR("AlignmentBase::swapTarget() Index out of range.", exception);
	swap(seqTemplateName[index1], seqTemplateName[index2]);
	swap(seqTemplate[index1], seqTemplate[index2]);
	swap(startAaTemplates[index1], startAaTemplates[index2]);
}


void
AlignmentBase::insertCharacter(unsigned int p, char c)
{
	if (p <= target.size())
		target.insert(p, 1, c);

	for (unsigned int i = 0; i < seqTemplate.size(); ++i)
		if (p <= seqTemplate[i].size())
			seqTemplate[i].insert(p, 1, c);
}


void
AlignmentBase::insertDash(unsigned int p)
{
	insertCharacter(p, '-');
}


void
AlignmentBase::deletePos(unsigned int p)
{
	PRECOND((p < target.size()), exception);
	target = deleteChar(target, p);

	for (unsigned int i = 0; i < seqTemplate.size(); ++i)
	{
		ASSERT((p < seqTemplate[i].size()), exception);
		seqTemplate[i] = deleteChar(seqTemplate[i], p);
	}
}


void
AlignmentBase::purgeTargetInsertions()
{
	unsigned int i = 0;

	while (i < target.size())
		if ((target[i] == '-') || (target[i] == 'X'))
		{
			deletePos(i);
			i = 0;
		}
		else
			++i;
}


void
AlignmentBase::cutTemplate(unsigned int index)
{
	// Removes all templates below index.
	if (index > seqTemplate.size())
		ERROR("Index out of scope.", exception);

	// Resize templates.
	seqTemplateName.resize(index);
	seqTemplate.resize(index);
	startAaTemplates.resize(index);

	// Cut empty positions.
	unsigned int i = 0;

	while (i < target.size())
		if (target[i] == '-')
		{
			bool gap = true;
			for (unsigned int j = 0; j < seqTemplate.size(); j++)
				gap = gap && (seqTemplate[j][i] == '-');
			if (gap)
				deletePos(i);
			else
				i++;
		}
		else
			i++;
}


/// The algorithm is so long, because 6 cases are distinguished:
/// - leading insertion for target A or target B
/// - normal insertion for target A or target B
/// - trailing insertion for target A or target B
/// If ignoreInsertion is set, insertions of target A are ignored.
void
AlignmentBase::addAlignment(const AlignmentBase &orig)
{
	unsigned int i = 0;
	unsigned int origPos = 0;
	unsigned int newTemplatePos = 0;
	unsigned int leftGapCount1 = 0;
	unsigned int leftGapCount2 = 0;
	unsigned int rightGapCount1 = 0;
	unsigned int rightGapCount2 = 0;
	unsigned int leftGapCount2Orig = leftGapCount2;
	unsigned int rightGapCount2Orig = rightGapCount2;
	int diff = 0;
	AlignmentBase other = orig; // make safety copy
	string pure1 = getPureSequence(target); // sequence without '-'
	string pure2 = getPureSequence(other.target); // sequence without '-'

	if (target == other.target)
		goto COMBINE_END; // nothing to do if sequences have same insertions

	// Prepare sequence for adding (provisorically add 'Z' for undefined aa).
	if (orig.startAaTarget > startAaTarget)
		for (int i = 0; i < (orig.startAaTarget - startAaTarget); ++i)
			other.insertCharacter(0, '~');
	else
		if (orig.startAaTarget < startAaTarget)
		{
			for (int i = 0; i < (orig.startAaTarget - startAaTarget); ++i)
				insertCharacter(0, '~');
			startAaTarget = orig.startAaTarget;
		}

	// Count still missing letters and add at end (getPureSequence also counts '-').
	diff = static_cast<int>(getPureSequence(target).size()) -
		static_cast<int>(getPureSequence(other.target).size());
	if (diff < 0) // still not same size: add "X" at end
		for (int i = 0; i < abs(diff); ++i)
			insertCharacter(target.size(), '~');
	else
		if (diff > 0) // still not same size: add "X" at end
			for (int i = 0; i < abs(diff); ++i)
				other.insertCharacter(other.target.size(), '~');

	#ifdef DEBUG_VERBOSE
	cout << "After preprocessing:\n"
	     << target << "\n"
	     << other.target << endl;
	#endif

	// ASSERT(getPureSequence(target).size() == getPureSequence(other.target).size(), exception);

	// Count leading gaps.
	while ((leftGapCount1 < target.size()) && (target[leftGapCount1] == '-'))
		++leftGapCount1;
	while ((leftGapCount2 < other.target.size()) && (other.target[leftGapCount2] == '-'))
		++leftGapCount2;
	leftGapCount2Orig = leftGapCount2;

	#ifdef DEBUG_VERBOSE
	cout << "Initial gap 1: " << leftGapCount1 << "\n"
	     << "Initial gap 2: " << leftGapCount2 << endl;
	#endif

	// Count trailing gaps.
	while ((rightGapCount1 < target.size()) && (target[target.size()-rightGapCount1-1] =='-'))
		++rightGapCount1;
	while ((rightGapCount2 < other.target.size()) && (other.target[other.target.size()-rightGapCount2-1] == '-'))
		++rightGapCount2;
	rightGapCount2Orig = rightGapCount2;

	#ifdef DEBUG_VERBOSE
	cout << "Trailing gap 1: " << rightGapCount1 << "\n"
	     << "Trailing gap 2: " << rightGapCount2 << endl;
	#endif

	if (leftGapCount1 > leftGapCount2)
	{
		for (unsigned int k = 0; k < (leftGapCount1 - leftGapCount2); ++k)
		{
			#ifdef DEBUG_VERBOSE
			cout << "Inserting gap into pos 0 of target 2: " << endl;
			#endif

			other.insertDash(0);
		}
		leftGapCount2 = leftGapCount1;
	}
	else
		if (leftGapCount1 < leftGapCount2)
		{
			for (unsigned int k = 0; k < (leftGapCount2 - leftGapCount1); ++k)
			{
				#ifdef DEBUG_VERBOSE
				cout << "Inserting gap into pos 0 of target 1: " << endl;
				#endif

				insertDash(0);
			}
			leftGapCount1 = leftGapCount2;
		}

	#ifdef DEBUG_VERBOSE
	cout << target << "\n"
	     << other.target << endl;
	#endif

	if (rightGapCount1 > rightGapCount2)
	{
		for (unsigned int k = 0; k < (rightGapCount1 - rightGapCount2); ++k)
		{
			#ifdef DEBUG_VERBOSE
			cout << "Inserting gap into pos " << target.size() << " of target 2: " << endl;
			#endif

			other.insertDash(other.target.size());
		}
		rightGapCount2 = rightGapCount1;
	}
	else
		if (rightGapCount1 < rightGapCount2)
		{
			for (unsigned int k = 0; k < (rightGapCount2 - rightGapCount1); ++k)
			{
				#ifdef DEBUG_VERBOSE
				cout << "Inserting gap into pos " << target.size() << " of target 1: " << endl;
				#endif

				insertDash(target.size());
			}
			rightGapCount1 = rightGapCount2;
		}

	#ifdef DEBUG_VERBOSE
	cout << target << "\n" << other.target << endl;
	#endif

	i = leftGapCount1; // now skip leading gaps
	while (i < (target.size() - rightGapCount1))
	{
		if (target[i] == '-')
		{
			origPos = getOrigPos(target, i); // position of first non-dash char
			newTemplatePos = getNewPos(other.target,origPos); // compute position on other sequence

			unsigned int l = i + 1; // count length of insertion
			for (l = i + 1; (l < target.size()) && (target[l] == '-'); ++l)
			{ }

			unsigned int l2 = l - i; // length of insertion
			unsigned int l3 = newTemplatePos;

			#ifdef DEBUG_VERBOSE
			cout << "Insertion of length " << l2 << " at position " << i << " found in target 1 " << "\n"
			     << "Corresponding position on target 2: " << l3 << endl;
			#endif

			for (l3 = newTemplatePos; (l3 < other.target.size()) && (other.target[l3] == '-'); ++l)
			{ }

			unsigned int l4 = l3 - newTemplatePos; // length of insertion of template

			#ifdef DEBUG_VERBOSE
			cout << "Found insertion of length " << l4 << " at position " << newTemplatePos << " of target 2" << endl;
			#endif

			if (l4 < l2) // insert this many dashes
			{
				unsigned int dl = l2 - l4;

				for (unsigned int j = 0; j < dl; ++j) // insertion of other alignment is dl residues too short
				{
					#ifdef DEBUG_VERBOSE
					cout << " inserting dash after position " << (newTemplatePos + 1) << " on target 2 " << endl;
					#endif

					other.insertDash(newTemplatePos+1);
				}

				i += dl;
			}
			else
				++i;
		}
		else
			++i;
	}

	i = leftGapCount2Orig; // start again from zero

	#ifdef DEBUG_VERBOSE
	cout << "Other direction:" << endl;
	#endif

	while (i < (orig.target.size() - rightGapCount2Orig)) // go through all '-' regions in target 2
	{
		if (orig.target[i] == '-')
		{
			#ifdef DEBUG_VERBOSE
			cout << "character " << i << " of target 2\n"
			     << " dash found." << endl;
			#endif

			origPos = getOrigPos(orig.target, i);
			newTemplatePos = getNewPos(target, origPos); // compute position on other sequence

			unsigned int l = i + 1; // count length of insertion
			for (l = i + 1; (l < orig.target.size()) && (orig.target[l] == '-'); ++l)
			{ }

			unsigned int l2 = l - i; // length of insertion
			unsigned int l3 = newTemplatePos; // count length of insertion at other alignment

			for (l3 = newTemplatePos; (l3 < target.size())&& (target[l3] == '-'); ++l)
			{ }

			unsigned int l4 = l3 - newTemplatePos; // length of insertion of template

			#ifdef DEBUG_VERBOSE
			cout << "Found insertion of length " << l2 << " at position " << i << " of target 2\n"
			     << "Found corresponding insertion of length " << l4 << " at position " << newTemplatePos << " of target 1" << endl;
			#endif

			if (l4 < l2) // insert this many dashes
			{
				unsigned int dl = l2 - l4;

				for (unsigned int j = 0; j < dl; ++j)// insertion of other alignement is dl residues too short
				{
				#ifdef DEBUG_VERBOSE
				cout << "Inserting dash into template 1 at position " << (newTemplatePos + 1) << endl;
				#endif

				insertDash(newTemplatePos+1);
				}

				i += dl;
			}
			else
				++i;
		}
		else
			++i;
	}

	COMBINE_END:

	#ifdef DEBUG_VERBOSE
	cout << "Result:\n"
	     << target << "\n"
	     << other.target << endl;
	#endif

	ASSERT((target.size() == other.target.size()), exception);

	// Append other templates.
	for (unsigned int i = 0; i < other.seqTemplate.size(); ++i)
	{
		seqTemplateName.push_back(other.seqTemplateName[i]);
		seqTemplate.push_back(other.seqTemplate[i]);
	}

	// Do clean up (replace '~' with '-').
	for (unsigned int i = 0; i < target.size(); ++i)
		if (target[i] == '~')
			target[i] = '-';

	for (unsigned int j = 0; j < seqTemplate.size(); ++j)
		for (unsigned int i = 0; i < seqTemplate[j].size(); ++i)
			if (seqTemplate[j][i] == '~')
				seqTemplate[j][i] = '-';
}


// HELPERS:

string
AlignmentBase::getPureSequence(const string &s)
{
	string result = "";

	for (unsigned int i = 0; i < s.size(); ++i)
		if (s[i] != '-')
			result.append(1, s[i]);

	return result;
}


/// If original index points to a dash, it returns the position of the first
/// non-dash left character, again taking '-' not into account.
unsigned int
AlignmentBase::getOrigPos(const string &s, unsigned int p)
{
	PRECOND((p < s.size()), exception);

	int i;

	for (i = static_cast<int>(p); i >= 0; --i) // move pointer to first non-dash character
		if (s[i] != '-')
			break;

	if (i < 0)
	{
//		ERROR("No defined original position.", exception);
		cout << "getOrigPos.i = " << i << endl;
		return 0; // left bound dashes
	}

	unsigned int counter = 0;

	for (int j = 0; j < i; ++j) // count the dashes before this position
		if (s[j] == '-')
			++counter;

	return i - counter; // position minus number of dashes counted
}


unsigned int
AlignmentBase::getNewPos(const string &s, unsigned int p)
{
	PRECOND((p < s.size()), exception);

	int count = -1; // count non-dash characters
	int pi = static_cast<int>(p);

	for (unsigned int i = 0; i < s.size(); ++i)
		if (s[i] != '-')
		{
			++count;
			if (count == pi)
				return i;
		}

	ERROR("Could not assign new index for sequence with insertions.", exception);

	return s.size(); // dummy
}


vector<string>
AlignmentBase::getTokens(const string &text)
{
	istringstream ist(text.c_str());
	char *charLine = new char[text.size()+1]; // size of string
	vector<string> v;
	string s;

	while (!ist.eof())
	{
		ist >> charLine;
		s = charLine; // assignment of c-strings to string
//		DUMP(s);
		if (s != "")
			v.push_back(s);
	}

	delete[] charLine;

	return v;
}


string
AlignmentBase::deleteChar(const string &s, unsigned int n)
{
	string result = s.substr(0, n);

	if (n < (s.size() - 1))
		result = result + s.substr(n+1, s.size());

	return result;
}

} // namespace
