#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

typedef unsigned long long uint64;

// 米勒罗宾素数测试算法(具有一定出错概率)
const uint64 Rabin_prime[5] = {2, 3, 5, 233, 331};

uint64 qmul(uint64 x, uint64 y, uint64 mod)
{
    return (x * y - (uint64)(x / (long double)mod * y + 1e-3) *mod + mod) % mod;// 乘法防止溢出， 如果p * p不爆uint64的话可以直接乘； O(1)乘法或者转化成二进制加法
}

uint64 qpow(uint64 a, uint64 n, uint64 mod)
{
    uint64 ret = 1;
    while (n)
    {
        if (n & 1)
        {
            ret = qmul(ret, a, mod);
        }
        a = qmul(a, a, mod);
        n >>= 1;
    }
    return ret;
}

bool Miller_Rabin(uint64 p)
{
    if (p < 2)
    {
        return false;
    }
    if (p != 2 && p % 2 == 0)
    {
        return false;
    }
    uint64 s = p - 1;
    while (! (s & 1))
    {
        s >>= 1;
    }
    for (uint64 i = 0; i < 5; ++i)
    {
        if (p == Rabin_prime[i])
        {
            return true;
        }
        uint64 t = s, m = qpow(Rabin_prime[i], s, p);
        while (t != p - 1 && m != 1 && m != p - 1)
        {
            m = qmul(m, m, p);
            t <<= 1;
        }
        if (m != p - 1 && !(t & 1))
        {
            return false;
        }
    }
    return true;
}

// 判断x是否为质数
bool judge_prime(uint64 x)
{
    if (x < 2)
    {
        return false;
    }
    for (uint64 cur_divnum = 2; cur_divnum <= sqrt((long double)x); ++cur_divnum)
    {
        if (x % cur_divnum == 0)
        {
            return false;
        }
    }
    return true;
}

// 计算产生从2-n之间所有的质数
vector<uint64> generate_prime(uint64 n)
{
    vector<uint64> prime;
    for (uint64 cur_num = 2; cur_num <= n; ++cur_num)
    {
        if (judge_prime(cur_num))
        // if (Miller_Rabin(cur_num))
        {
            prime.push_back(cur_num);
        }
    }
    return prime;
}

// 通过二分查找找到n在prime向量中的第一个不小于n的数的下标
// prime为从小到大排列的质数向量，无重复数字
uint64 find_end_index(uint64 n, const vector<uint64>& prime)
{
    uint64 start_index = 0;
    uint64 end_index = prime.size() - 1;
    uint64 mid_index = (start_index + end_index)/2;
    while (start_index < end_index)
    {
        mid_index = (start_index + end_index)/2;
        if (n == prime[mid_index])
        {
            return mid_index;
        }
        else if (n > prime[mid_index])
        {
            if (mid_index + 1 < end_index)
            {
                if (n > prime[mid_index + 1])
                {
                    start_index = mid_index + 1;
                }
                else
                {
                    return mid_index + 1;
                }
            }
            else
            {
                return end_index;
            }
        }
        else if (n < prime[mid_index])
        {
            if (mid_index >= start_index + 1)// +1在右边，若放到左边变为-1，当mid_index为0时，左侧会溢出
            {
                if (n <= prime[mid_index - 1])
                {
                    end_index = mid_index - 1;
                }
                else
                {
                    return mid_index;
                }
            }
            else
            {
                return mid_index;
            }
        }
    }
    return end_index;// 若n>prime[prime.size() - 1],则返回prime.size() - 1
}

// 进行因式分解
vector<uint64> factorize(uint64 n, const vector<uint64>& prime)
{
    vector<uint64> factorize_number;
    if (judge_prime(n))// 如果n为质数，无须分解，直接返回
    // if (Miller_Rabin(n))// 如果n为质数，无须分解，直接返回
    {
        factorize_number.push_back(n);
        return factorize_number;
    }
    uint64 cur_num = n;
    uint64 cur_index = 0;
    uint64 end_index = find_end_index(sqrt((long double)cur_num), prime);
    while (cur_index <= end_index)
    {
        if (cur_num % prime[cur_index] == 0)
        {
            factorize_number.push_back(prime[cur_index]);
            cur_num /= prime[cur_index];
            end_index = find_end_index(sqrt((long double)cur_num), prime);
        }
        else
        {
            cur_index++;
        }
    }
    if (cur_num != 1)// 如果最后剩的质数不是1，则把最后一个质数也添加进去
    {
        factorize_number.push_back(cur_num);
    }
    return factorize_number;
}

// 显示vector<uint64>的内容
void print_vector(const vector<uint64>& v)
{
    for (vector<uint64>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        cout << *it << " ";
    }
    cout << endl;
}

// 显示分解质因数的结果
void print_factorize_number(uint64 n, const vector<uint64>& v)
{
    if (v.size() == 1)// 该数为质数
    {
        cout << n << " = 1 * "<< v[0];
    }
    else
    {
        vector<uint64>::const_iterator it = v.begin();
        cout << n << " = " << *it;
        for (++it; it != v.end(); ++it)
        {
            cout << " * " << *it;
        }
    }
    cout << endl;
}

void test1()
{
    uint64 test_number = 7140229933;
    vector<uint64> prime = generate_prime(sqrt((long double)test_number));
    // print_vector(prime);
    print_factorize_number(test_number, factorize((long double)test_number, prime));
}

void test2()
{
    uint64 start_number = 6541367000;
    uint64 end_number = 6541367999;
    vector<uint64> prime = generate_prime(sqrt((long double)end_number));
    // print_vector(prime);
    for (uint64 cur_number = start_number; cur_number <= end_number; ++cur_number)
    {
        print_factorize_number(cur_number, factorize((long double)cur_number, prime));
    }
}

void test3()
{
    uint64 start_number = 6541367000;
    uint64 end_number = 6541367999;
    vector<uint64> prime = generate_prime(sqrt((long double)end_number));
    // print_vector(prime);
    for (uint64 cur_number = start_number; cur_number <= end_number; ++cur_number)
    {
        vector<uint64> factorize_result = factorize((long double)cur_number, prime);
        if (factorize_result.size() == 2)
        {
            print_factorize_number(cur_number, factorize_result);
        }
    }
}

int main(int argc, char** argv)
{
    test1();
    // test2();
    // test3();
    system("pause");
    return 0;
}
