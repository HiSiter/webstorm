var biRadixBits = 16;
var biRadix = 1 << 16; // = 2^16 = 65536


//取模*
function exp_mod(a, n, b) {
    var t;
    if (n == 0) return 1 % b;
    if (n == 1) return a % b;
    t = exp_mod(a, n / 2, b);
    t = t * t % b;
    if ((n & 1) == 1) t = t * a % b;
    return t;
}

//判断是否是质数
function isprime(num) {
    if (num == 1) return 0;
    if (num == 2) return 1;

    var i = 2;
    for (i; i <= Math.sqrt(num); i++) {
        if (num % i == 0) {
            return 0;
            break;
        }
    }
    return 1;
}

//判断两个数是否互质
function gcd(x, y) {
    var max, min, temp;
    var x = parseInt(x);
    var y = parseInt(y);
    max = x > y ? x : y;
    min = x < y ? x : y;
    while (max % min) {

        temp = max % min;
        max = min;
        min = temp;
    }
    return min;
}

//大数相加
function bigNumAdd(a, b) {
    var m = a.split('').reverse();
    var n = b.split('').reverse();
    var ret = [];
    var s = 0;

    for (var i = 0; i < a.length || i < b.length; i++) {
        var t = (m[i] | 0) + (n[i] | 0) + s;

        ret.push(t % 10);
        s = (t / 10) | 0;
    }
    if (s) {
        ret.push(s);
    }
    return ret.reverse().join('');
}

//异或运算
function xor(a, b) {
    var res = a ^ b;
    return res;
}


//大数相乘
function bigNumMulti(a, b) {
    var p = a.match(/\d{1,4}/g).reverse();
    var q = b.match(/\d{1,4}/g).reverse();
    var f1 = 0;
    var result = "0";

    for (var i = 0; i < p.length; i++) {
        var f2 = 0;
        for (var j = 0; j < q.length; j++) {
            var t = (p[i] | 0) * (q[j] | 0);
            t += new Array(f1 + f2 + 1).join("0");
            result = bigNumAdd(result, t);
            f2 += q[j].length;
        }
        f1 += p[i].length;
    }
    return result;
}

//求次方
function bigNumPow(a, b) {
    var ret = "1";
    for (var i = 0; i < b; i++) {
        ret = bigNumMulti(ret, a.toString());
    }

    return ret;
}


//乘法
function cheng(a, b) {
    var zong = [String(a), String(b)];
    var fen = [];
    zong = getMax(zong[0], zong[1]);

    zong[0] = zong[0].split('');
    zong[1] = zong[1].split('');
    //获取b的长度,处理乘法分配率的乘法
    for (var j = (zong[1].length - 1); j >= 0; j--) {
        var next = 0;
        var fentemp = [];
        var jialing = '';
        //获取a的长度处理乘法
        for (var i = (zong[0].length - 1); i >= 0; i--) {
            var ji = Number(zong[0][i]) * Number(zong[1][j]) + next;
            fentemp.unshift(ji % 10);
            next = Math.floor(ji / 10);
            if (i == 0 && !(next == 0)) {
                fentemp.unshift(next);
            }
        }
        //后面添加0
        jialing = new Array((zong[1].length - (j + 1)) + 1).join('0');
        fentemp.push(jialing);
        fen[j] = fentemp.join('');
    }
}

//找逆
function inverse(X, L) {
    var X = parseInt(X);
    var L = parseInt(L);
    var temp = 1;
    var temp1 = 0;
    while (temp <= L[0]) {
        //if (temp.multiply(X).mod(L).equals(one)) {
        temp1 = cheng(temp + "", X + "");
        temp1 = temp1 % L[0];
        //temp1 = temp;
        //temp = (temp * X) % L;
        if (temp1 == 1) {
            return temp;
        }
        temp++;
    }
    return temp;
}


function biSubtract(x, y) {
    var result;
    if (x.isNeg != y.isNeg) {
        y.isNeg = !y.isNeg;
        result = biAdd(x, y);
        y.isNeg = !y.isNeg;
    } else {
        result = new BigInt();
        var n, c;
        c = 0;
        for (var i = 0; i < x.digits.length; ++i) {
            n = x.digits[i] - y.digits[i] + c;
            result.digits[i] = n & 0xffff;
            // Stupid non-conforming modulus operation.
            if (result.digits[i] < 0) result.digits[i] += biRadix;
            c = 0 - Number(n < 0);
        }
        // Fix up the negative sign, if any.
        if (c == -1) {
            c = 0;
            for (var i = 0; i < x.digits.length; ++i) {
                n = 0 - result.digits[i] + c;
                result.digits[i] = n & 0xffff;
                // Stupid non-conforming modulus operation.
                if (result.digits[i] < 0) result.digits[i] += biRadix;
                c = 0 - Number(n < 0);
            }
            // Result is opposite sign of arguments.
            result.isNeg = !x.isNeg;
        } else {
            // Result is same sign.
            result.isNeg = x.isNeg;
        }
    }
    return result;
}

function biMultiplyMod(x, y, m) {
    return biModulo(biMultiply(x, y), m);
}

//遍历求逆
function getInverseOfGk(k, p) {
    var res;
    for (var i = 1; i < 999999999; i++) {
        if ((i * parseInt(k) - 1) % p == 0) {
            res = i;
            return res;
        }
    }

}

//计算结果
function showResult() {
    var result = "";
    var p = form.p.value;
    var x = form.x.value;
    var g = form.g.value;
    var M = form.M.value;
    var k = form.k.value;


    var y;
    var s;
    var newM;
    var r;
    // var JustTest = parseInt(getInverseOfGk(29, 29));

    // alert("test=" + JustTest);


    //判断p是否为质数
    if (isprime(p) != 1) {
        result += "p 必须是一个质数";
        return;

    }
    //判断x与p-1是否互质
    else if (gcd(x, p - 1) == 1) {
        // var temp = bigNumPow(g, x);

        //y = g^x mod p
        y = parseInt(exp_mod(g, x, p));

        var temp = parseInt(M) % (p - 1);
        // alert(temp);

        //求s ； s=（y+M）^(M mod (p-1)) mod p
        s = exp_mod((parseInt(y) + parseInt(M)), parseInt(temp), parseInt(p));
        //
        var ni;
        var gkmodp = exp_mod(g, k, p);
        var ni = getInverseOfGk(gkmodp, p);
        //
        var t1 = parseInt(M) % p;
        var t2 = parseInt(s) % p;
        //求r
        r = (t1 * t2 * ni) % p;
        var XORrs = xor(r, s)
        //x^-1：x在p-1的逆元
        var x_1mop = getInverseOfGk(x, p - 1);
        //s^-1：s在p的逆元
        var s_1mop = getInverseOfGk(s, p);

        //求s+t
        var sAddt;
        var e1, e2;
        e1 = x_1mop % (p - 1);
        e2 = (k % (p - 1) - XORrs % (p - 1) + (p - 1)) % (p - 1);


        sAddt = (e1 * e2) % (p - 1);
        //求M'
        newM = parseInt(M % p);
        result += "1.\t y = " + y + "\n";

        result += "2.\t s = " + s + "\n";

        result += "3\t g^k * g^(-k) = 1 (mod p) " + "\t\t";

        result += gkmodp + " * " + ni + " = 1 (mod " + p + ")\n";

        result += "4.\t r = " + r + "\n";

        result += "5.\t x * x^(-1) = 1 (mod (p-1))" + "\t\t";

        result += x + " * " + x_1mop + " = 1 (mod (" + p + "-1))" + "\n";

        result += "6.\t Xor (r,s) =" + XORrs + "\n";

        result += "7.\t s + t = " + sAddt + "\n";

        result += "8.\t s * s^(-1) = 1 (mod p)" + "\t\t";

        result += s + " * " + s_1mop + " = 1 (mod " + p + ")" + "\n";

        result += "9.\t M' = " + newM + "\n";

    }

    form.final_result.value = (result);
    result = ""
}

