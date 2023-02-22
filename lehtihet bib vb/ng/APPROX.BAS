Attribute VB_Name = "APPROX"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << APPROX >>
'  Requis   :           << EVAL >>
'  Externes :
'
'  Descrip. : Approximations polynomiales
'
Option Explicit
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'          Interpolation polynomiale simple
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Function Lagrange(ByVal k1 As Integer, _
                  ByVal k2 As Integer, _
                  ByRef x() As Double, _
                  ByRef y() As Double, _
                  ByVal xo As Double) As Double
    
    Dim i As Integer, j As Integer
    Dim xi As Double, xj As Double
    Dim Pol As Double, Lag As Double
    
    If (k2 < k1) Then Stop ' Données non valides
    Pol = 0#
    For i = k1 To k2
        xi = x(i)
        If (xo = xi) Then Lagrange = y(i): Exit Function
        Lag = 1#
        For j = k1 To k2
            If i <> j Then
                xj = x(j)
                If xi = xj Then Stop ' Points non distincts
                Lag = Lag * ((xo - xj) / (xi - xj))
            End If
        Next j
        Pol = Pol + y(i) * Lag
    Next i
    Lagrange = Pol
End Function
Function Hermite(ByVal k1 As Integer, _
                 ByVal k2 As Integer, _
                 ByRef x() As Double, _
                 ByRef y() As Double, _
                 ByRef d() As Double, _
                 ByVal xo As Double) As Double
    
    Dim i As Integer, j As Integer
    Dim xi As Double, xj As Double
    Dim Pol As Double, Lag As Double
    Dim s As Double, h As Double
    
    If (k2 < k1) Then Stop
    Pol = 0#
    For i = k1 To k2
        xi = x(i)
        If (xo = xi) Then Hermite = y(i): Exit Function
        Lag = 1#: s = 0#
        For j = k1 To k2
            If (i <> j) Then
                xj = x(j)
                If xi = xj Then Stop ' Points non distincts
                h = 1 / (xi - xj)
                s = s + h
                Lag = Lag * ((xo - xj) * h)
            End If
        Next j
        h = xo - xi
        Pol = Pol + Lag * Lag * ((1# - 2# * s * h) * y(i) + h * d(i))
    Next i
    Hermite = Pol
End Function
Function Neville(ByVal k1 As Integer, _
                 ByVal k2 As Integer, _
                 ByRef x() As Double, _
                 ByRef y() As Double, _
                 ByVal xo As Double) As Double
    
    Dim k As Integer, M As Integer
    Dim N As Integer, ko As Integer, kk As Integer
    Dim w() As Double, h As Double, U As Double
    
    ko = k1 - 1
    N = k2 - ko
    If (N < 1) Then Stop ' Données non valides
    If (N = 1) Then Neville = y(k1): Exit Function
    ReDim w(N)
    For k = 1 To N
        w(k) = y(k + ko)
    Next k
    For M = 1 To N - 1
        For k = 1 To N - M
            kk = k + ko
            U = x(kk)
            h = U - x(M + kk)
            If (h = 0) Then Stop ' Points non distincts
            w(k) = w(k) + (w(k) - w(k + 1)) * ((xo - U) / h)
        Next k
    Next M
    Neville = w(1)
End Function
Function Newton(ByVal k1 As Integer, _
                ByVal k2 As Integer, _
                ByRef x() As Double, _
                ByRef c() As Double, _
                ByVal xo As Double) As Double
    
    Dim k As Integer, ko As Integer, N As Integer
    Dim Pol As Double
    
    N = k2 - k1 + 1
    If (N < 1) Then Stop ' Données non valides
    If N = 1 Then Newton = c(1): Exit Function
    Pol = c(N)
    ko = k2 - N
    For k = N - 1 To 1 Step -1
       Pol = (xo - x(k + ko)) * Pol + c(k)
    Next k
    Newton = Pol
End Function
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'        Interpolation polynomiale par morceaux
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Function Lineaire(ByVal N As Integer, _
                  ByRef x() As Double, _
                  ByRef y() As Double, _
                  ByVal xo As Double, _
                  ByRef yo As Double) As Double

    Dim k1 As Integer, k2 As Integer
    Dim h As Double, v As Double, w As Double

    Call COEFF.Morceau(N, x(), xo, k1, k2, h, v, w)
    Lineaire = w * y(k1) + v * y(k2)
End Function
Function Cubique(ByVal N As Integer, _
                 ByRef x() As Double, _
                 ByRef y() As Double, _
                 ByRef d() As Double, _
                 ByVal xo As Double, _
                 ByRef yo As Double) As Double

    Dim k1 As Integer, k2 As Integer
    Dim h As Double, v As Double, w As Double

    Call COEFF.Morceau(N, x(), xo, k1, k2, h, v, w, u1, u2)
    Cubique = w * w * (y(k1) * (1# + 2# * v) + u1 * d(k1)) _
            + v * v * (y(k2) * (1# + 2# * w) - u2 * d(k2))
End Function
Function Spline(ByVal N As Integer, _
                ByRef x() As Double, _
                ByRef y() As Double, _
                ByRef s() As Double, _
                ByVal xo As Double) As Double
                
    Dim k1 As Integer, k2 As Integer
    Dim h As Double, v As Double, w As Double

    Call COEFF.Morceau(N, x(), xo, k1, k2, h, v, w)
    Spline = w * y(k1) + v * y(k2) + _
            (w * (w * w - 1#) * s(k1) + _
             v * (v * v - 1#) * s(k2)) * (h * h / 6)
End Function
