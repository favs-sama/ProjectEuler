import math 

def multiples_of_three_and_five():
    sum = 0
    for i in range(1000):
        if i % 3 == 0 or i % 5 == 0:
            sum = sum + i 
    print(sum)
    
def even_fibonacci_numbers():
    prev = 1 
    curr = 2
    sum = curr 
    while curr < 4000000:
        t = curr 
        curr = curr + prev 
        prev = t 
        if curr % 2 == 0:
            sum = sum + curr
    print(sum)

# simple algorithm to extract prime factors   
def largest_prime_factor():
    number = 600851475143
    largest = 0
    sqrt = int(math.sqrt(number))
    for i in range(2, sqrt):
        if number % i == 0:
            largest = max(largest, i) 
            while number % i == 0:
                number = number / i
    largest = max(largest, number)
    print(largest)
    
# start from 999 to get to the largest palindrome
# faster 
def largest_palindrome():
    largest = 0
    for i in range(999, 99, -1):
        for j in range(i-1, 99, -1):
            number = str(i * j)
            if number == number[::-1]:
                if largest > i*j:
                    break
                largest = max(largest, i*j)
    print(largest)
    
#helper
def generate_prime_factors(n):
    prime_factors = []
    sqrt = int(math.sqrt(n)) + 1
    for i in range(2, sqrt):
        if n % i == 0:
            prime_factors.append(i)
            n = n / i
            while n % i == 0:
                prime_factors.append(i)
                n = n / i
    if n != 1:
        prime_factors.append(int(n))
    return prime_factors

# this can be solved via dynamic programming is the solution for n
# can be calculated by checking the difference in the set of prime
# factors between solution for n - 1 and number n, the difference 
# is then multiplied with the answer for n - 1
# dp[3] => 6 => [2, 3]
# dp[4] => ? => [2, 2]
# [2, 3] - [2, 2] = [2]
# dp[4] = [2, 3, 2] => 12 
# dp[4] = 12 => [2, 2, 3]
def smallest_multiple(n):
    prime_factors = [0, 1, 2]
    dp = [0, 1, 2]
    for i in range(3, n + 1):
        if dp[i - 1] % i == 0:
            dp.append(dp[i - 1])
        else:
            # check if the number can be formed from the 
            # existing prime factors
            diff = list(set(prime_factors).intersection(generate_prime_factors(i)))
            if len(diff) < 1:
                dp.append(dp[i - 1] * i)
                prime_factors.append(i)
            else:
                mult = 1
                for j in diff:
                    mult = mult * j
                    prime_factors.append(j)
                dp.append(dp[i - 1] * mult)
    return dp[n]

def sum_square_difference():
    squares = 0
    sums = 0 
    for i in range(1, 101):
        squares = squares + (i * i)
        sums = sums + i 
    return sums ** 2 - squares
    
    
def ten_thousand_one_prime(place):
    i = 2
    curr = 4
    while i < place + 1:
        isPrime = True 
        num = int(math.sqrt(curr)) + 1
        for p in range(2, num):
            if curr % p == 0:
                isPrime = False
                break
            
        if isPrime:
            i = i + 1
            if i == place:
                print(curr)
        curr = curr + 1

# add window for 13 numbers and get product and update max per iteration
def largest_product_in_a_series():
    series="7316717653133062491922511967442657474235534919493496983520312774506326239578318016984801869478851843858615607891129494954595017379583319528532088055111254069874715852386305071569329096329522744304355766896648950445244523161731856403098711121722383113622298934233803081353362766142828064444866452387493035890729629049156044077239071381051585930796086670172427121883998797908792274921901699720888093776657273330010533678812202354218097512545405947522435258490771167055601360483958644670632441572215539753697817977846174064955149290862569321978468622482839722413756570560574902614079729686524145351004748216637048440319989000889524345065854122758866688116427171479924442928230863465674813919123162824586178664583591245665294765456828489128831426076900422421902267105562632111110937054421750694165896040807198403850962455444362981230987879927244284909188845801561660979191338754992005240636899125607176060588611646710940507754100225698315520005593572972571636269561882670428252483600823257530420752963450"
    largest_product=0
    for i in range(987):
        subset=series[i:i+13]
        product=1
        for s in subset:
            product = product * int(s)
        largest_product = max(product, largest_product)
    print(largest_product)
    
def special_pythagorean_triplet():
    for i in range(1, 1000):
        for j in range(i, 1000):
            k = 1000 - i - j 
            if i**2 + j**2 == k**2:
                return i * j * k

# implement Sieve of Eratosthenes to save prime numbers        
def summation_of_primes(n):
    primes = [True for i in range(n + 1)]
    p = 2 
    while p ** 2 <= n:
        if primes[p] == True:
            for i in range(p**2, n + 1, p):
                primes[i] = False
        p = p + 1
    
    sums = 0
    for i in range(2, len(primes)):
        if primes[i] == True:
            sums = sums + i 
    return sums
    
def sieve_of_eratosthenes(n):
    # WIP - optimize
    primes = [False for i in range(int((n - 1) / 2) + 1)]
    crosslimit = int(math.sqrt(n + 1) - 1 / 2)
    for i in range(crosslimit):
        if primes[i] == False:
            for j in range(2 * i * (i + 1), int((n - 1) / 2 + 1), 2 * i + 1):
                primes[i] = True
    
    sums = 0 
    for i in range(2, len(primes)):
        if primes[i] == False:
            sums = sums + (2 * i + 1)
    return sums

#implement checking for all directions
def largest_product_in_a_grid():
    grid = []
    grid.append([8, 2, 22, 97, 38, 15, 0, 40, 0, 75, 4, 5, 7, 78, 52, 12, 50, 77, 91, 8])
    grid.append([49, 49, 99, 40, 17, 81, 18, 57, 60, 87, 17, 40, 98, 43, 69, 48, 4, 56, 62, 0])
    grid.append([81, 49, 31, 73, 55, 79, 14, 29, 93, 71, 40, 67, 53, 88, 30, 3, 49, 13, 36, 65])
    grid.append([52, 70, 95, 23, 4, 60, 11, 42, 69, 24, 68, 56, 1, 32, 56, 71, 37, 2, 36, 91])
    grid.append([22, 31, 16, 71, 51, 67, 63, 89, 41, 92, 36, 54, 22, 40, 40, 28, 66, 33, 13, 80])
    
    grid.append([24, 47, 32, 60, 99, 3, 45, 2, 44, 75, 33, 53, 78, 36, 84, 20, 35, 17, 12, 50])
    grid.append([32, 98, 81, 28, 64, 23, 67, 10, 26, 38, 40, 67, 59, 54, 70, 66, 18, 38, 64, 70])
    grid.append([67, 26, 20, 68, 2, 62, 12, 20, 95, 63, 94, 39, 63, 8, 40, 91, 66, 49, 94, 21])
    grid.append([24, 55, 58, 5, 66, 73, 99, 26, 97, 17, 78, 78, 96, 83, 14, 88, 34, 89, 63, 72])
    grid.append([21, 36, 23, 9, 75, 0, 76, 44, 20, 45, 35, 14, 0, 61, 33, 97, 34, 31, 33, 95])
    
    grid.append([78, 17, 53, 28, 22, 75, 31, 67, 15, 94, 3, 80, 4, 62, 16, 14, 9, 53, 56, 92])
    grid.append([16, 39, 5, 42, 96, 35, 31, 47, 55, 58, 88, 24, 0, 17, 54, 24, 36, 29, 85, 57])
    grid.append([86, 56, 0, 48, 35, 71, 89, 7, 5, 44, 44, 37, 44, 60, 21, 58, 51, 54, 17, 58])
    grid.append([19, 80, 81, 68, 5, 94, 47, 69, 28, 73, 92, 13, 86, 52, 17, 77, 4, 89, 55, 40])
    grid.append([4, 52, 8, 83, 97, 35, 99, 16, 7, 97, 57, 32, 16, 26, 26, 79, 33, 27, 98, 66])
    
    grid.append([88, 36, 68, 87, 57, 62, 20, 72, 3, 46, 33, 67, 46, 55, 12, 32, 63, 93, 53, 69])
    grid.append([4, 42, 16, 73, 38, 25, 39, 11, 24, 94, 72, 18, 8, 46, 29, 32, 40, 62, 76, 36])
    grid.append([20, 69, 36, 41, 72, 30, 23, 88, 34, 62, 99, 69, 82, 67, 59, 85, 74, 4, 36, 16])
    grid.append([20, 73, 35, 29, 78, 31, 90, 1, 74, 31, 49, 71, 48, 86, 81, 16, 23, 57, 5, 54])
    grid.append([1, 70, 54, 71, 83, 51, 54, 69, 16, 92, 33, 48, 61, 43, 52, 1, 89, 19, 67, 48])
    
    largest_product = 0
    
    #horizontal
    for row in grid:
        for i in range(17):
            product = row[i] * row[i + 1] * row[i + 2] * row[i + 3]
            largest_product = max(largest_product, product)
            
    #vertical
    for i in range(20):
        for j in range(17):
            product = grid[j][i] * grid[j + 1][i] * grid[j + 2][i] * grid[j + 3][i]
            largest_product = max(product, largest_product)

    #diagonal left-right
    #iterate through columns
    largest_diagonal = 0    

    upper_bound = 20 
    for i in range(20):
        if upper_bound < 4:
            break
        k = i
        for j in range(upper_bound - 3):
            product = grid[j][k] * grid[j + 1][k + 1] * grid[j + 2][k + 2] * grid[j + 3][k + 3]
            largest_product = max(largest_product, product)
            k = k + 1
        upper_bound = upper_bound - 1
        
    #print(largest_diagonal)
        
    #iterate through rows
    upper_bound = 19
    for i in range(1, 20):
        if upper_bound < 4:
            break
        k = i
        for j in range(upper_bound - 3):
            product = grid[k][j] + grid[k + 1][j + 1] + grid[k + 2][j + 2] + grid[k + 3][j + 3]
            largest_product = max(product, largest_product)
            k = k + 1
        upper_bound = upper_bound - 1
    
    #diagonal right-left
    upper_bound = 19
    for col in range(19, -1, -1):
        k = col 
        for j in range(upper_bound - 2):
            product = grid[j][k] * grid[j + 1][k - 1] * grid[j + 2][k - 2] * grid[j + 3][k - 3]
            largest_product = max(product, largest_product)
            k = k - 1
        upper_bound = upper_bound - 1
        
    upper_bound = 18
    for row in range(1, 20):
        if upper_bound < 3:
            break
        k = row 
        for j in range(19, 19 - upper_bound + 2, -1):
            product = grid[k][j] * grid[k + 1][j - 1] * grid[k + 2][j - 2] * grid[k + 3][j - 3]
            largest_product = max(product, largest_product)
            k = k + 1
        upper_bound = upper_bound - 1

    print(largest_product)
    
def arithmetic_series(n):
    return int(n * (1 + n) / 2)

def highly_divisible_triangular_number():
    nth = 3 
    triangular_number = 0 
    divisor = 0 
    while divisor <= 500:
        triangular_number = arithmetic_series(nth)
        divisor = 2
        for i in (2, int(math.sqrt(triangular_number)) + 1):
            if triangular_number % i == 0:
                divisor = divisor + 2 
        nth = nth + 1
  
#can be optimized with DP      
def longest_collatz_sequence():
    starting_number = 0
    max_sequence_length = 0 
    for i in range(2, 1000000):
        sequence_length = 1
        j = i 
        while j > 1:
            if j % 2 == 0:
                j = j / 2
            else:
                j = 3 * j + 1
            sequence_length = sequence_length + 1
        if sequence_length > max_sequence_length:
            starting_number = i 
            max_sequence_length = sequence_length
            
    print(starting_number)

dp = {}
dp[1] = 1
dp[2] = 2 

# recursive approach with dynamic programming since 
# a lot of sequences are repeated at some point so it 
# is better to save it for use at some point
def count_chain(n):
    if n in dp:
        return dp[n]
        
    if n % 2 == 0:
        dp[n] = 1 + count_chain(n / 2)
    else:
        dp[n] = 2 + count_chain(int((3 * n + 1) / 2))
    
    return dp[n]
    

def longest_collatz_sequence_dp():
    longest_chain = 0
    number = 0 
    for i in range(3, 1000000):
        count = count_chain(i)
        if count > longest_chain:
            longest_chain = count 
            number = i 
            
    print(number)
 
# combinatorics this is simply n! / (k! * (n - k)!)  
def lattice_paths():
    factorial = [0, 1, 2]
    
    for i in range(3, 41):
        factorial.append(i * factorial[i - 1])

    return factorial[40] / (factorial[20] ** 2)
    
    
    
def power_digit_sum():
    power = str(1 << 1000)
    pds = 0 
    for i in range(len(power)):
        pds = pds + int(power[i])
        
    print(pds)
    
def number_letter_counts():
    #counting numbers
    lengths = {0 : 0, 1 : 3, 2 : 3, 3 : 5, 4 : 4, 5 : 4, 6 : 3, 7 : 5, 8 : 5, 9 : 4 }
    
    #teens
    lengths[11] = 6 
    lengths[12] = 6 
    lengths[13] = 8 
    lengths[14] = 8 
    lengths[15] = 7
    lengths[16] = 7 
    lengths[17] = 9 
    lengths[18] = 8 
    lengths[19] = 8
    
    #under hundreds edges
    lengths[10] = 3 #ten
    lengths[20] = 6 #
    lengths[30] = 6
    lengths[40] = 5  
    lengths[50] = 5 
    lengths[60] = 5 
    lengths[70] = 7 
    lengths[80] = 6 
    lengths[90] = 6 
    
    count = 0
    
    for i in range(1, 20):
        count = count + lengths[i]
        
    #under hundreds save to dictionary to retrieve in thousands
    for i in range(20, 100):
        if i < 30 and i % 20 <= 9:
            lengths[i] = lengths[20] + lengths[i % 20]
        elif i < 40 and i % 30 <= 9:
            lengths[i] = lengths[30] + lengths[i % 30]
        elif i < 50 and i % 40 <= 9:
            lengths[i] = lengths[40] + lengths[i % 40]
        elif i < 60 and i % 50 <= 9:
            lengths[i] = lengths[50] + lengths[i % 50]
        elif i < 70 and i % 60 <= 9:
            lengths[i] = lengths[60] + lengths[i % 60]
        elif i < 80 and i % 70 <= 9:
            lengths[i] = lengths[70] + lengths[i % 70]
        elif i < 90 and i % 80 <= 9:
            lengths[i] = lengths[80] + lengths[i % 80]
        elif i < 100 and i % 90 <= 9:
            lengths[i] = lengths[90] + lengths[i % 90]
        count = count + lengths[i]
        
    count = count + 10
    for i in range(101, 1000):
        count = count + 10 + lengths[int(i / 100)] + lengths[i % (int(i / 100) * 100)]
            
    count = count + 11 #for a thousand
    count = count - 24 #eliminate and in 100s/200s
        
    print(count)

#update each level by getting maximum sum for each subtree
#almost similar to minimum path sum when offset
def maximum_path_sum():
    triangle = []
    triangle.append([75])
    triangle.append([95, 64])
    triangle.append([17, 47, 82])
    triangle.append([18, 35, 87, 10])
    triangle.append([20, 4, 82, 47, 65])
    triangle.append([19, 1, 23, 75, 3, 34])
    triangle.append([88, 2, 77, 73, 7, 63, 67])
    triangle.append([99, 65, 4, 28, 6, 16, 70, 92])
    triangle.append([41, 41, 26, 56, 83, 40, 80, 70, 33])
    triangle.append([41, 48, 72, 33, 47, 32, 37, 16, 94, 29])
    triangle.append([53, 71, 44, 65, 25, 43, 91, 52, 97, 51, 14])
    triangle.append([70, 11, 33, 28, 77, 73, 17, 78, 39, 68, 17, 57])
    triangle.append([91, 71, 52, 38, 17, 14, 91, 43, 58, 50, 27, 29, 48])
    triangle.append([63, 66, 4, 68, 89, 53, 67, 30, 73, 16, 69, 87, 40, 31])
    triangle.append([4, 62, 98, 27, 23, 9, 70, 98, 73, 93, 38, 53, 60, 4, 23])
    
    for i in range(13, -1, -1):
        for j in range(len(triangle[i])):
            triangle[i][j] = triangle[i][j] +  max(triangle[i + 1][j], triangle[i + 1][j + 1])
            
    return triangle[0][0]

#the idea here is in a perfect world where every month starts 
#in sunday, there would be 28 days each month. Since this is 
#not a perfect world each month offsets the start day depending
#on the remainder when the days of the month are divided by 7
def counting_sundays():
    months = {1 : 31, 2 : 28, 3 : 31, 4 : 30, 5 : 31, 6 : 30, 7 : 31, 8 : 31, 9 : 30, 10 : 31, 11 : 30, 12 : 31}
    start_day = 2 
    sundays = 0
    for i in range(1901, 2001):
        if i % 4 == 0:
            months[2] = 29 
        else :
            months[2] = 28 
        
        for j in range(1, 13):
            if start_day == 0:
                sundays = sundays + 1 
            start_day = (start_day + (months[j] % 7)) % 7 
            
    print(sundays)
    
def factorial_digit_sum():
    factorial = 1
    for i in range(1, 101):
       factorial = factorial * i 
       
    f = str(factorial)
    
    digit_sum = 0 
    for i in range(len(f)):
        digit_sum = digit_sum + int(f[i])
        
    print(digit_sum)
    
def amicable_numbers():
    #sieve of eratosthenes to save prime numbers
    primes = [True for i in range(10001)]
    p = 2 
    while p ** 2 <= 10000:
        if primes[p] == True:
            for i in range(p**2, 10001, p):
                primes[i] = False
        p = p + 1

    d = {}
    #save in a dictionary the sum of divisors
    for i in range(1, 10000):
        if (primes[i] == False):
            d[i] = 0
            for j in range(2, int(math.sqrt(i) + 1)):
                if i % j == 0:
                    d[i] = d[i] + j + int(i / j)
            d[i] = d[i] + 1


    sum_amicable_numbers = 0
    #check the amicable numbers
    for key in d: 
        if d[key] in d and d[d[key]] == key and d[key] != d[d[key]]:
            sum_amicable_numbers = sum_amicable_numbers + key
            
    print(sum_amicable_numbers)
    
def IsAbundant(n):
    if n == 1:
        return False
    i = 2 
    total = 1
    while (i < math.ceil(math.sqrt(n))):
        if n % i == 0:
            total = total + i + (n // i)
        i = i + 1 
        
    if i ** 2 == n:
        total = total + i 
        
    if total > n:
        return True 
    return False
    
    #sieve of eratosthenes to save prime numbers
    primes = [True for i in range(10001)]
    p = 2 
    while p ** 2 <= 10000:
        if primes[p] == True:
            for i in range(p**2, 10001, p):
                primes[i] = False
        p = p + 1
        
def non_abundant_sums():
    abundant_numbers = []
    for i in range(12, 28124):
        if IsAbundant(i):
            abundant_numbers.append(i)
            
    two_abundant = [False for i in range(28124)]
    
    for i in range(len(abundant_numbers)):
        for j in range(0, len(abundant_numbers), 1):
            if abundant_numbers[i] + abundant_numbers[j] <= 28123:
                two_abundant[abundant_numbers[i] + abundant_numbers[j]] = True 
            
    sums = 0     
    for i in range(28124):
        if not two_abundant[i]:
            sums = sums + i 
            
    print(sums)
    
#compute for factoriadic number as this
#indicates position of numbers in the nth
#permutation (actually n -1th)
def decimal_to_factoriadic(n):
    i = 1 
    result = ""
    while n > 0:
        result = str(n % i) + result
        n = int(n / i)
        i = i + 1
        
    return result
    
def lexicographic_permutations(place):
    numbers = [i for i in range(10)]
    factoriadic = decimal_to_factoriadic(place)
    answer = ""
    for i in factoriadic:
        answer = answer + str(numbers[int(i)])
        numbers.remove(numbers[int(i)])
    return answer
    
#time and memory efficient fibonacci
def nth_fibonacci_with_n_digits(n):
    first = 1 
    second = 1 
    i = 2
    while len(str(second)) < n:
        first, second = second, first + second
        i = i + 1 
    return i

#input copy pasted in console :3
def names_scores():
    names = input()
    names_list = names.split(',')
    
    names_list_modified = []
    for n in names_list:
        names_list_modified.append(n.strip('"'))
    
    names_list_modified.sort()
    
    letters = {'A' : 1, 'B' : 2, 'C' : 3, 'D' : 4, 'E' : 5, 'F' : 6, 'G' : 7, 'H' : 8, 'I' : 9, 'J' : 10, 'K' : 11, 'L' : 12, 'M' : 13, 'N' : 14, 'O' : 15, 'P' : 16, 'Q' : 17, 'R' : 18, 'S' : 19, 'T' : 20, 'U' : 21, 'V' : 22, 'W' : 23, 'X' : 24, 'Y' : 25, 'Z' : 26 }
    
    score = 0;
    for i in range(len(names_list_modified)):
        local_sum = 0 
        name = names_list_modified[i]
        for j in name:
            local_sum = local_sum + letters[j]
        score = score + ((i + 1) * local_sum)
    
    #answer should be 871198282
    print(score)
        
    
print(nth_fibonacci_with_n_digits(1000))
    


        
    
















