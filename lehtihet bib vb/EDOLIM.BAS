Attribute VB_Name = "EDOLIM"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << EDOLIM >>
'  Requis   :           << SYSLIN >>
'  Externes :  Fonc(nf,x)
'
'  Descrip. : Equations diff. Ordinaires
'             Problèmes de conditions aux limites
'
Option Explicit
Function DFC(ByVal nf As Integer, _
             ByVal nfR As Integer, _
             ByVal nfD As Integer, _
             ByVal a As Double, _
             ByVal b As Double, _
             ByVal Va As Double, _
             ByVal Vb As Double, _
             ByVal N As Integer, _
             ByRef x() As Double, _
             ByRef y() As Double) As Integer

    Dim k As Integer, m As Integer
    Dim u As Double, w As Double
    Dim h As Double, h2 As Double, hh As Double
    Dim Bas() As Double, Haut() As Double, Diag() As Double
    
    If N < 2 Then Stop
    ReDim x(0 To N), y(0 To N)
    m = N - 1
    ReDim Bas(1 To m), Haut(1 To m), Diag(1 To m)
    h = (b - a) / N
    hh = h * h: h2 = 0.5 * h
    u = a
    For k = 1 To m
        u = u + h
        x(k) = u
        y(k) = hh * fonc(nf, u)
        w = h2 * fonc(nfD, u)
        Bas(k) = 1# - w
        Haut(k) = 1# + w
        Diag(k) = -2# + hh * fonc(nfR, u)
    Next k
    x(0) = a: y(0) = Va
    x(N) = b: y(N) = Vb
    y(1) = y(1) - Bas(1) * y(0)
    y(m) = y(m) - Haut(m) * y(N)
    Debug.Print Diag(1), Haut(1), 0, y(1)
    Debug.Print Bas(2), Diag(2), Haut(2), y(2)
    Debug.Print 0, Bas(3), Diag(3), y(3)
        
    DFC = SYSLIN.Thomas(m, Bas(), Diag(), Haut(), y())
End Function

