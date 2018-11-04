/**
 * This is an implementation of Huffman coding.
 *
 * The core algorithm is taken from the CLR book (Introduction of Algorithms),
 * Chapter 16.3, and directly used to implement the 'build_tree()' routine.
 *
 * After the tree is built, a code table that maps a character to a binary
 * code is built from the tree, and used for encoding text. Decoding is done
 * by traversing the Huffman tree, as prescribed by the algorithm.
 *
 * Binary codes are represented by std::vector<bool>, which is a specialized
 * vector that optimizes space.
 */

#include <vector>
#include <queue>
#include <map>
#include <algorithm>
#include <string>
#include <cassert>
#include <stdexcept>
#include <iostream>

using namespace std;

// A Huffman Tree Node
//哈夫曼树的结点
struct HuffmanTree {
	//字符
	char c; // character in an alphabet
	//字符对应出现的次数
	int cfreq; // frequency of c.
	struct HuffmanTree *left;
	struct HuffmanTree *right;
	HuffmanTree(char c, int cfreq, struct HuffmanTree *left = NULL,
		struct HuffmanTree *right = NULL) :
		c(c), cfreq(cfreq), left(left), right(right) {
	}
	~HuffmanTree() {
		delete left, delete right;
	}
	//比较两个结点大小的方法
	// Compare two tree nodes
	class Compare {
	public:
		bool operator()(HuffmanTree *a, HuffmanTree *b) {
			return a->cfreq > b->cfreq;
		}
	};
};

/**
 * Builds a Huffman Tree from an input of alphabet C, where C is a vector
 * of (character, frequency) pairs.
 */
//创建一棵哈夫曼树
//通过一个字符与出现次数的pair的向量，生成哈夫曼树的结点,从而生成树
HuffmanTree *build_tree(vector< pair<char, unsigned> > &alph) {
	// First build a min-heap
	// Build leaf nodes first
	//优先队列--最小堆--用来存储树的结点
	//字符出现的次数越小，优先级越高
	priority_queue<HuffmanTree *, vector<HuffmanTree *>, HuffmanTree::Compare > alph_heap;
	//遍历向量，生成树的结点
	//并将结点push入最小堆中
	for (vector< pair<char, unsigned> >::iterator it = alph.begin();
		it != alph.end(); ++it) {
		//这里生成的是叶子结点，所以left和right都是NULL
		HuffmanTree *leaf = new HuffmanTree(it->first, it->second);
		alph_heap.push(leaf);
	}

	// HuffmanTree algorithm: Merge two lowest weight leaf nodes until
	// only one node is left (root).
	//每次生成的新树根
	HuffmanTree *root = NULL;
	//将最小堆中最小的两个进行合并，生成新的父结点
	while (alph_heap.size() > 1) {
		HuffmanTree *l, *r;
		l = alph_heap.top();
		alph_heap.pop();
		r = alph_heap.top();
		alph_heap.pop();
		root = new HuffmanTree(0, l->cfreq + r->cfreq, l, r);
		//将新生出来的结点push入最小堆中
		alph_heap.push(root);
	}
	//返回根节点
	return root;
}

/**
 * Prints the tree nodes in breadth-first order
 */
//广度优先打印哈夫曼树
void print_tree(HuffmanTree *t) {
	//双端队列--对于队列两端的操作效率较高
	//下面的int表示这个结点位于树的第几个层次
	deque< pair<HuffmanTree *, int> > q;
	//根节点位于第0层
	q.push_back(make_pair(t, 0));
	int curlevel = -1;
	while (!q.empty()) {
		//广度优先遍历--每次都是从队头出来
		//如果存在左右孩子，就将左右孩子push入队尾
		HuffmanTree *parent = q.front().first;
		//得到该结点在第几层
		int level = q.front().second;
		q.pop_front();
		//只有当层数变化的时候，才打印层数
		if (curlevel != level) {
			curlevel = level;
			cout << "Level " << curlevel << endl;
		}

		//打印这个结点的字符以及出现的次数
		cout << parent->cfreq << " " << parent->c << endl;
		if (parent->left)
			q.push_back(make_pair(parent->left, level + 1));
		if (parent->right)
			q.push_back(make_pair(parent->right, level + 1));
	}
}

//从哈夫曼树中创建一个映射表
//将一个字符映射成一个二进制数
//这里使用一个bool的向量表示二进制数
//01使用的是bool类型表示的
typedef vector<bool> code_t;
typedef map<char, code_t> codetable;
/**
 * Makes a lookup table (std::map) of (c -> code) from a HuffmanTree, where
 * code is an unsigned long representing the binary code.
 */
 //返回的是一个map里面包含了字符及其对应的二进制编号
map<char, code_t> build_lookup_table(HuffmanTree *htree) {
	//最终的映射表
	codetable lookup;
	//双端队列--存放的是结点以及它对应的二进制编号
	deque< pair<HuffmanTree *, code_t> > q;
	//根节点--创建一个空的bool向量
	q.push_back(make_pair(htree, code_t()));
	//广度优先--创建每个叶子节点的编号
	while (!q.empty()) {
		HuffmanTree *node, *lc, *rc;
		code_t code;
		node = q.front().first;
		code = q.front().second;
		q.pop_front();
		lc = node->left;
		rc = node->right;
		//哈夫曼树总是慢慢的--没有节点或者是有两个结点
		if (lc) {
			// HuffmanTree is always full (either no children or two children)
			// Left child is appended a 0 and right child a 1.
			//得到上一次的编码
			code_t code_cp(code);
			//在上一次的基础上加0或者是加1
			q.push_back(make_pair(lc, (code.push_back(0), code)));
			q.push_back(make_pair(rc, (code_cp.push_back(1), code_cp)));
		}
		//叶子结点--需要保存字符类型
		else {
			// Leaf node: contains the character
			lookup.insert(make_pair(node->c, code));
			cout << "(" << node->c << ", ";
			//遍历bool类型的向量，打印出这个字符对应的编码
			for (unsigned i = 0; i < code.size(); i++) {
				cout << code[i];
			}
			cout << ")" << endl;
		}
	}

	return lookup;
}

/**
 * Encodes an input string. returns a byte vector.
 */
 //输入一个string,以及密码本，输出它的哈夫曼编码
code_t encode(string input, codetable &lookup) {
	//bool类型的数组--二进制编码
	code_t result;

	//对string进行遍历
	for (string::iterator it = input.begin(); it != input.end(); ++it) {
		//通过map中*it这个字符，查找它对应的bool类型的向量
		code_t b = lookup[*it];
		//将得到的二进制编码添加到最终的bool类型的向量后面
		result.insert(result.end(), b.begin(), b.end());
	}
	//将编码后的结果向量返回
	return result;
}

/**
 * Look up the next valid code in @biter using @htree and returns the
 * resulting string. Note the iterator @biter is advanced by the actual
 * length of the next valid code, which varies.
 */
 //通过输入的二进制编码是0还是1决定往左走还是往右走，
 //走到叶子节点的时候将叶子节点的字符返回
char code_lookup(code_t::iterator &biter, const code_t::iterator &biter_end,
	const HuffmanTree *htree) {
	const HuffmanTree *node = htree;

	while (true) {
		//哈夫曼树总是满满的--有两个结点或者没有结点
		//node是一个叶子节点
		if (!node->left) {
			// Huffman tree is full: always contains both children or none.
			//跳出循环--将叶子节点对应的字符返回
			break;
		}
		//二进制编码匹配结束
		if (biter == biter_end) {
			throw std::out_of_range("No more bits");
		}
		//还有二进制编码
		//并且为1--向右走
		if (*biter) {
			node = node->right;
		}
		//二进制编码为0--向左走
		else {
			node = node->left;
		}
		++biter;
	}

	return node->c;
}

/**
 * Decodes a compressed string represented by a bit vector (vector<char>)
 * @compressed, using a HuffmanTree @htree.
 * Returns the original string.
 */
 //解码过程
 //输入二进制编码，以及哈夫曼树,输出一个字符串
string decode(code_t &compressed, const HuffmanTree *htree) {
	string result;

	code_t::iterator biter = compressed.begin();
	//注意下面传的是引用，所以iterator的位置真的会改变的
	while (true) {
		try {
			result += code_lookup(biter, compressed.end(), htree);
		}
		catch (const std::out_of_range &oor) {
			// Iterator exhausted.
			break;
		}
	}

	return result;
}

/**
 * Tests
 */
 // Make frequency table from a string.
//测试
//通过string构造一个字符与出现次数映射关系的向量
vector< pair<char, unsigned> > make_freq_table(string inp) {
	map<char, unsigned> cfmap;
	vector< pair<char, unsigned> >cfvec;

	//遍历整个string
	//如果这个字符在map中没有出现过就创建一个pair
	//如果出现过就进行++
	for (unsigned i = 0; i < inp.size(); i++) {
		//找到了结尾还是没有找到，就证明这个字符在map中没有出现过
		if (cfmap.find(inp[i]) == cfmap.end()) {
			cfmap.insert(make_pair(inp[i], 1));
		}
		cfmap[inp[i]] += 1;
	}

	//遍历map将里面的pair放入向量中
	for (map<char, unsigned>::iterator it = cfmap.begin();
		it != cfmap.end(); ++it) {
		cfvec.push_back(make_pair(it->first, it->second));
	}

	return cfvec;
}
//下面不理解
string bitvec_to_string(code_t &bitvec) {
	string result;
	size_t nbits;

	nbits = bitvec.size() & 7;

	// Write the number of "hanging bits" at the first byte
	result += static_cast<char>(nbits); // at most 7

	char byte = 0;
	for (unsigned i = 0; i < bitvec.size(); i++) {
		unsigned boff = i & 7;
		byte |= bitvec[i] << boff;
		if (boff == 7) {
			// Write a byte
			result += byte;
			byte = 0;
		}
	}
	if (nbits) {
		result += byte;
	}

	return result;
}

code_t string_to_bitvec(string packed) {
	code_t result;

	assert(packed.size());
	if (packed.size() == 1) {
		assert(packed[0] == 0);
		return result;
	}
	unsigned nbits = packed[0];
	for (string::iterator it = packed.begin() + 1; it != packed.end(); ++it) {
		for (unsigned i = 0; i < 8; i++) {
			result.push_back((*it >> i) & 1);
		}
	}
	// fix the last byte
	if (nbits) {
		for (unsigned i = 0; i < (8 - nbits); i++) {
			result.pop_back();
		}
	}

	return result;
}

#include <cstdio>
void hexdump(const unsigned char *bytes, int nbytes) {
	int i, j;

	for (i = 0; i < nbytes; i += 16) {
		printf("%06x: ", i);
		for (j = 0; j < 16; j++) {
			if (i + j < nbytes) {
				printf("%02x ", bytes[i + j]);
			}
			else {
				printf("   ");
			}
		}
		printf(" ");
		for (j = 0; j < 16; j++)
			if (i + j < nbytes)
				printf("%c", isprint(bytes[i + j]) ? bytes[i + j] : '.');
		printf("\n");
	}
}

string tests[] = {
	"aDDzzz",
	"Declaration of Independence",
	"We hold these truths to be self-evident, that all men are created equal, that they are endowed by their Creator with certain unalienable Rights, that among these are Life, Liberty and the pursuit of Happiness.--That to secure these rights, Governments are instituted among Men, deriving their just powers from the consent of the governed, --That whenever any Form of Government becomes destructive of these ends, it is the Right of the People to alter or to abolish it, and to institute new Government, laying its foundation on such principles and organizing its powers in such form, as to them shall seem most likely to effect their Safety and Happiness.",
	"int main() {for (unsigned i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {string s = tests[i]; cout << \"\\n>>>>>>> test \" << i << \" >>>>>>>\" << endl; vector< pair<char, unsigned> > cfvec = make_freq_table(s); HuffmanTree *htree = build_tree(cfvec); //print_tree(htree); codetable ctbl = build_lookup_table(htree); code_t t = encode(s, ctbl); cout << \"original:\" << endl << s << endl; cout << \"encoded (compression ratio: \" << (t.size() + 7) / 8 << \"/\" << s.size() << \" or \" << ((float)(t.size() / 8) / s.size()) << \"):\" << endl; string packed = bitvec_to_string(t); hexdump((unsigned char *)(packed.c_str()), packed.size()); // Decode code_t t1 = string_to_bitvec(packed); assert(std::equal(t.begin(), t.end(), t1.begin())); string s1 = decode(t1, htree); cout << \"decoded:\\n\" << s1 << endl; assert(s1 == s); delete htree;}}"
};

int main() {
	for (unsigned i = 0; i < sizeof(tests) / sizeof(tests[0]); i++) {
		string s = tests[i];
		cout << "\n>>>>>>> test " << i << " >>>>>>>" << endl;
		//通过字符串生成一个字符于其对应次数的向量
		vector< pair<char, unsigned> > cfvec = make_freq_table(s);
		//通过字符与其次数的映射向量创造出一棵哈夫曼树
		HuffmanTree *htree = build_tree(cfvec);
		//print_tree(htree);
		//通过哈夫曼树创建一个密码表（字符与其对应的哈夫曼编码的映射）
		codetable ctbl = build_lookup_table(htree);
		//通过stirng以及对应的密码表得到整个string的哈夫曼编码
		code_t t = encode(s, ctbl);
		cout << "original:" << endl << s << endl;
		//计算压缩的比率
		cout << "encoded (compression ratio: "
			<< (t.size() + 7) / 8 << "/" << s.size() << " or "
			<< ((float)(t.size() / 8) / s.size()) << "):" << endl;
		//string packed = bitvec_to_string(t);
		//cout << "packed:" << packed << endl;
		//cout << "t:";
		//for (code_t::iterator it = t.begin(); it != t.end(); it++) {
		//	cout << *it;
		//}
		//cout << endl;
		//hexdump((unsigned char *)(packed.c_str()), packed.size());

		//// Decode解码
		//code_t t1 = string_to_bitvec(packed);
		//for (code_t::iterator it = t1.begin(); it != t1.end(); it++) {
		//	cout << *it;
		//}
		//cout << endl;
		//assert(std::equal(t.begin(), t.end(), t1.begin()));
		//通过整个字符串的二进制编码以及哈夫曼树进行解码
		//string s1 = decode(t1, htree);
		string s1 = decode(t, htree);
		cout << "decoded:" << s1 << endl;
		assert(s1 == s);
		delete htree;
	}
	system("pause");
	return 0;
}