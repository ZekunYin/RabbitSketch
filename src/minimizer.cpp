#include "minimizer.h"
#include <unordered_map>
using namespace std;

//hash64 -from minimap2
uint64_t hash64(uint64_t key, uint64_t mask){
	key = (~key + (key << 21)) & mask;//key = (key << 21) - key - 1;
	key = key ^ key>>24;
	key = ((key + (key << 3)) + (key << 8)) & mask;//key * 265
	key = key ^ key>>14;
	key = ((key + (key << 2)) + (key << 4)) & mask;//key * 21
	key = key ^ key>>28;
	key = (key + (key << 31)) & mask;
	return key;
}

bool findElement(vector <uint64_t> arr, uint64_t element){
	bool result =false;
	for(int i = 0; i < arr.size(); i++){
		if(arr[i] == element){
			result = true;
			break;
		}
	}
	return result;
}



//vector <uint64_t> findMinimizers(int k, int w, string s)
void findMinimizers(int k, int w, string s, vector <uint64_t> &minimizerSketch)
{

	uint8_t seq_nt4_table[256] = {
		0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
		4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
	};

	if(w < 0 || w > 256){
		std::cerr << "w must be: 0 < w < 257" << std::endl;
		exit(1);
	}
	if(k < 0 || k > 31){
		std::cerr << " k size must be: 0 < k < 32" << std::endl;	
		exit(1);
	}

	int len = s.length();
	//cout << "the len is: " << len << endl;
	if(len < 1){
		std::cerr << "sequence length must be > 0 " << std::endl;
		exit(1);
	}
	if(len < (w + k - 1)){
		std::cerr << "sequence length must be >= w + k -1" << std::endl;
		exit(1);
	}
//	vector <uint64_t> minimizerSketch;

	uint64_t kmers[2] = {0, 0};
	int kmerSpan = 0;

	//cout << "the k is: " << k << endl;
	uint64_t bitmask = (uint64_t)(((uint64_t)1 << (2 * k)) - 1);
	uint64_t bitshift = (uint64_t)(2 * (k - 1));
//	printf("the bitmask is: %lx \n", bitmask);
//	printf("the bitshift is: %lx \n", bitshift);


	//q = queue.newQueue();
	vector <Pair> q;
	std::unordered_map<uint64_t, double> minimizerMap;

	for(int i = 0; i < len; i++){
		int windowIndex = i - w + 1;
		//cout << "the windowIndex is: " << windowIndex << endl; //done@xxm

		uint8_t c = seq_nt4_table[s[i]]; // c should be uint8_t
		//printf("%x\n",c);
		
		//printf("the c is: %d \n", c); //done@xxm
		if(c > 3){
			//TODO: handle these, by skiping the base and starting w again? @hulk

		}
		if(windowIndex + 1 < k){
			kmerSpan = windowIndex + 1;
		}
		else{
			kmerSpan = k;
		}

		//get the forward k-mer
		kmers[0] = (kmers[0] << 2 | (uint64_t)c) & bitmask;

		//get the reverse k-mer
		kmers[1] = (kmers[1] >> 2) | ((uint64_t)3 ^ (uint64_t)c) << bitshift;

		//don't try for minimizers until a full k-mer has been collected
		if(i < k - 1){
			//cout << "i < k - 1" << endl;
			continue;
		}

		//skip symmetric k-mers as we don't know the strand
		if(kmers[0] == kmers[1]){
			//cout << "kmers0 == kmers1" << endl;
			continue;
		}

		//set the canonical k-mer
		int strand = 0;
		
		if(kmers[0] > kmers[1]){
			strand = 1;
		}

		Pair currentKmer;

		uint64_t x = hash64(kmers[strand], bitmask) << 8 | (uint64_t)kmerSpan;
		//printf("%lx\n",kmers[strand]);
		
		//cout << "the i is: " << i << " and the hash64 x is: " << x << endl;
		//printf("%lx\n",x);
		currentKmer.X = x;
		currentKmer.Y = i;

		//if(q.size() != 0){
			//if minimizers are in the q from the previous window, remove them.
			while(q.size() != 0){
				if(q.front().Y > (i - w)){
					break;
				}
				q.erase(q.begin());
			}

			// hashed k-mers less than equal to the currentKmer are not required,so remove them from the back of the q.
			while(q.size() != 0){
				if(q.back().X < currentKmer.X){
					break;
				}
				q.pop_back();
			}
		//}

		q.push_back(currentKmer);
		//no problem before this.@xxm
	//	for(int i_ = 0; i_ < q.size(); i_++){
	//		printf("%lx\t",q[i_].X);
	//	}
	//	printf("\n");
		
		//printf("%lx\n",currentKmer.X);

		if(windowIndex >= 0){
		//	if(minimizerSketch.size() == 0){
		//		minimizerSketch.push_back(q.front().X);
		//		continue;
		//	}

		//	if(findElement(minimizerSketch, q.front().X)){
		//		continue;
		//	}

		//	minimizerSketch.push_back(q.front().X);
			minimizerMap.insert({q.front().X, 1.0});
		//	if(!findElement(minimizerSketch, q.front().X)){
		//		minimizerSketch.push_back(q.front().X);
		//		//printf("%lx\n",q.front().X);
		//	}

		}

	}//end of sequence.
	for(auto& x:minimizerMap){
		minimizerSketch.push_back(x.first);
	}

//	return minimizerSketch;
}

