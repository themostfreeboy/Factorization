#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

typedef unsigned long long uint64;

// �����ޱ����������㷨(����һ���������)
const uint64 Rabin_prime[5] = {2, 3, 5, 233, 331};

uint64 qmul(uint64 x, uint64 y, uint64 mod)
{
    return (x * y - (uint64)(x / (long double)mod * y + 1e-3) *mod + mod) % mod;// �˷���ֹ����� ���p * p����uint64�Ļ�����ֱ�ӳˣ� O(1)�˷�����ת���ɶ����Ƽӷ�
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

// �ж�x�Ƿ�Ϊ����
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

// ���������2-n֮�����е�����
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

// ͨ�����ֲ����ҵ�n��prime�����еĵ�һ����С��n�������±�
// primeΪ��С�������е��������������ظ�����
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
            if (mid_index >= start_index + 1)// +1���ұߣ����ŵ���߱�Ϊ-1����mid_indexΪ0ʱ���������
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
    return end_index;// ��n>prime[prime.size() - 1],�򷵻�prime.size() - 1
}

// ������ʽ�ֽ�
vector<uint64> factorize(uint64 n, const vector<uint64>& prime)
{
    vector<uint64> factorize_number;
    if (judge_prime(n))// ���nΪ����������ֽ⣬ֱ�ӷ���
    // if (Miller_Rabin(n))// ���nΪ����������ֽ⣬ֱ�ӷ���
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
    if (cur_num != 1)// ������ʣ����������1��������һ������Ҳ��ӽ�ȥ
    {
        factorize_number.push_back(cur_num);
    }
    return factorize_number;
}

// ��ʾvector<uint64>������
void print_vector(const vector<uint64>& v)
{
    for (vector<uint64>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        cout << *it << " ";
    }
    cout << endl;
}

// ��ʾ�ֽ��������Ľ��
void print_factorize_number(uint64 n, const vector<uint64>& v)
{
    if (v.size() == 1)// ����Ϊ����
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
