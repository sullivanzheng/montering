#include <iostream>

class counter{
//A counter with "lap" function as in a stopwatch.
private:
	unsigned long num;
	unsigned long lastreset;
public:
	counter(void){num=lastreset=0;}
	void clear(void){num=lastreset=0;}
	//inline void count(void){num++;}
	inline void operator++ (int){num++;}
	inline void lap(void){lastreset=num;}
    inline unsigned long getNumber(void){return num-lastreset;}
	inline unsigned long operator () (void){return num-lastreset;}
	inline unsigned long getTotCounts(void){return num;}
};

template <typename TypeName> class statqueue{
//A simple stat tool to give running min, max, mean, var, std, std of mean.
private:
    long itemNum;
    TypeName squareSum;
    TypeName sum,
        maxItem,minItem,lastItem;
public:
    statqueue(){
        itemNum=0;
        squareSum=sum=0;
        maxItem=minItem=0;
    }
    void clear(){
        itemNum=0;
        squareSum=sum=0;
        lastItem=maxItem=minItem=0;
	}
    void push(TypeName r_data){
        itemNum++;
		lastItem=r_data;
        sum+=r_data;
        squareSum+=(r_data*r_data);
        if (itemNum==1){
            maxItem=minItem=r_data;
        }
        else{
            maxItem=r_data>maxItem?r_data:maxItem;
            minItem=r_data<minItem?r_data:minItem;
        }
    }
    TypeName getMean(){
        return sum/itemNum;
    }
    TypeName getVar(){
        if (itemNum<=1) {
			std::cout<<"statequeue:Could not return variance when itemNum<=1"<<std::endl;
			getchar();exit(EXIT_FAILURE);
        }
        return (squareSum-(sum*sum)/itemNum)/(itemNum-1);
    }
    TypeName getStdev(){
        return sqrt(getVar());
    }
    TypeName getStdevOfMean(){
        return sqrt(getVar())
            /sqrt(double(itemNum));
    }
    TypeName getmaxItem(){
        if (itemNum<1) {
            std::cout<<"statequeue:Could not return min:max since it is empty"<<endl;
            getchar();exit(EXIT_FAILURE);
        }
        return maxItem;
    }
    TypeName getminItem(){
        if (itemNum<1) {
            std::cout<<"statequeue:Could not return min:max since it is empty"<<endl;
            getchar();exit(EXIT_FAILURE);
        }
        return minItem;
    }
	TypeName getlastItem(){
		if (itemNum==0){
			throw ("the queue is emply");
		}
		return lastItem;
	}
	TypeName getitemNum(){
		return itemNum;
	}
	TypeName getsquareSum(){
		return squareSum;
	}
	TypeName getsum(){
		return sum;
	}
};