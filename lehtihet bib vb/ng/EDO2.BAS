Attribute VB_Name = "EDO2"
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << EDO2 >>
'  Requis   :           << SYSLIN >>
'  Externes :            FoncPQR(nf,x,p,q,r)
'
'  Descrip. : Equations diff. Ordinaires
'
Option Explicit
Function DFC(ByVal nf As Integer, _
             ByVal a As Double, _
             ByVal b As Double, _
             ByVal ya As Double, _
             ByVal yb As Double, _
             ByVal N As Integer, _
             ByRef x() As Double, _
             ByRef y() As Double) As Integer
    
    ' Résolution de l'équation différentielle
    '
    '      y" = P(x)y'(x) + Q(x)y(x) + R(x)
    '
    ' par la méthode des différences finies
    
    Dim k As Integer, m As Integer
    Dim t As Double, w As Double
    Dim h As Double, h2 As Double, hh As Double
    Dim p As Double, q As Double, r As Double
    Dim L() As Double, U() As Double, D() As Double
    
    If N < 3 Then Stop
    m = N - 1
    ReDim x(0 To N), y(0 To N)
    ReDim L(1 To m), U(1 To m), D(1 To m)
    h = (b - a) / N
    hh = h * h:    h2 = 0.5 * h
    w = a
    For k = 1 To m
        w = w + h
        Call foncPQR(nf, w, p, q, r)
        x(k) = w
        y(k) = hh * r
        D(k) = -2# - hh * q
        t = h2 * p
        L(k) = 1# + t
        U(k) = 1# - t
    Next k
    x(0) = a: y(0) = ya
    x(N) = b: y(N) = yb
    y(1) = y(1) - L(1) * ya
    y(m) = y(m) - U(m) * yb
    DFC = SYSLIN.Thomas2(m, D(), L(), U(), y())
End Function



