package main

import (
	"math/big"
)

func eratosthenes(x int) (primes []int) {
	table := make([]bool, x+1)
	table[0] = true
	table[1] = true

	for i := 0; i <= x; i++ {
		if table[i] {
			continue
		}
		if i != 2 {
			primes = append(primes, i)
		}

		for j := i * i; j <= x; j += i {
			table[j] = true
		}
	}
	return
}

func legendre(a, p *big.Int) (b int) {
	temp := new(big.Int)
	temp.Div(temp.Sub(p, big.NewInt(1)), big.NewInt(2))
	temp.Exp(a, temp, p)
	b = int(temp.Int64())
	if b == 1 || b == 0 {
		return 1
	} else {
		return -1
	}
}

/*
func main() {
	a := big.NewInt(6469)
	b := big.NewInt(7919)
	fmt.Print(modSqrt(a, b))
}
*/
