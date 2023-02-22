Attribute VB_Name = "EDO"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << EDO >>
'  Requis   :
'  Externes :  EDO1_F(nf,t,y)
'             *EDO2_F(nf,t,y,z,Py,Pz)
'             *EDOS_F(nf,t,y(),p())
'
'  Descrip. : Equations diff. Ordinaires
'
Option Explicit
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'      Equations différentielles du premier ordre
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Sub EDO1_Euler(ByVal nf As Integer, _
               ByVal a As Double, _
               ByVal b As Double, _
               ByVal Ya As Double, _
               ByVal N As Integer, _
               ByVal m As Integer, _
               ByRef T() As Double, _
               ByRef y() As Double)

    Dim i As Integer, j As Integer
    Dim h As Double
    Dim u As Double, v As Double
    
    If N < 1 Or m < 1 Then Stop
    ReDim T(0 To N), y(0 To N)
    h = ((b - a) / N) / m
    u = a: v = Ya: i = 0
    Do
        T(i) = u: y(i) = v
        i = i + 1: If i > N Then Exit Sub
        For j = 1 To m
            v = v + h * EDO1_F(nf, u, v)
            u = u + h
        Next j
    Loop
End Sub
Sub EDO1_Heun(ByVal nf As Integer, _
              ByVal a As Double, _
              ByVal b As Double, _
              ByVal Ya As Double, _
              ByVal N As Integer, _
              ByVal m As Integer, _
              ByRef T() As Double, _
              ByRef y() As Double)

    Dim i As Integer, j As Integer
    Dim h As Double, h2 As Double
    Dim u As Double, v As Double
    Dim P0 As Double
    
    If N < 1 Or m < 1 Then Stop
    ReDim T(0 To N), y(0 To N)
    h = ((b - a) / N) / m
    h2 = 0.5 * h
    u = a: v = Ya: i = 0
    Do
        T(i) = u: y(i) = v
        i = i + 1: If i > N Then Exit Sub
        For j = 1 To m
            P0 = EDO1_F(nf, u, v)
            u = u + h
            v = v + h2 * (P0 + EDO1_F(nf, u, v + h * P0))
        Next j
    Loop
End Sub
Sub EDO1_PM(ByVal nf As Integer, _
            ByVal a As Double, _
            ByVal b As Double, _
            ByVal Ya As Double, _
            ByVal N As Integer, _
            ByVal m As Integer, _
            ByRef T() As Double, _
            ByRef y() As Double)

    Dim i As Integer, j As Integer
    Dim h As Double, h2 As Double
    Dim u As Double, v As Double
    Dim P0 As Double
    
    If N < 1 Or m < 1 Then Stop
    ReDim T(0 To N), y(0 To N)
    h = ((b - a) / N) / m
    h2 = 0.5 * h
    u = a: v = Ya: i = 0
    Do
        T(i) = u: y(i) = v
        i = i + 1:  If i > N Then Exit Sub
        For j = 1 To m
            P0 = EDO1_F(nf, u, v)
            u = u + h2
            v = v + h * EDO1_F(nf, u, v + h2 * P0)
            u = u + h2
        Next j
    Loop
End Sub
Sub EDO1_RK3(ByVal nf As Integer, _
             ByVal a As Double, _
             ByVal b As Double, _
             ByVal Ya As Double, _
             ByVal N As Integer, _
             ByVal m As Integer, _
             ByRef T() As Double, _
             ByRef y() As Double)

    Dim i As Integer, j As Integer
    Dim h As Double, h3 As Double, h4 As Double, h23 As Double
    Dim u As Double, v As Double
    Dim P0 As Double, Pm1 As Double, Pm2 As Double
    
    If N < 1 Or m < 1 Then Stop
    ReDim T(0 To N), y(0 To N)
    h = ((b - a) / N) / m
    h4 = h / 4#: h3 = h / 3#: h23 = h3 + h3
    u = a: v = Ya: i = 0
    Do
        T(i) = u:  y(i) = v
        i = i + 1:  If i > N Then Exit Sub
        For j = 1 To m
            P0 = EDO1_F(nf, u, v)
            Pm1 = EDO1_F(nf, u + h3, v + h3 * P0)
            Pm2 = EDO1_F(nf, u + h23, v + h23 * Pm1)
            v = v + h4 * (P0 + 3# * Pm2)
            u = u + h
        Next j
    Loop
End Sub
Sub EDO1_RK4(ByVal nf As Integer, _
             ByVal a As Double, _
             ByVal b As Double, _
             ByVal Ya As Double, _
             ByVal N As Integer, _
             ByVal m As Integer, _
             ByRef T() As Double, _
             ByRef y() As Double)

    Dim i As Integer, j As Integer
    Dim h As Double, h2 As Double, h6 As Double
    Dim u As Double, v As Double
    Dim P0 As Double, Pm1 As Double, Pm2 As Double, p1 As Double
    
    If N < 1 Or m < 1 Then Stop
    ReDim T(0 To N), y(0 To N)
    h = ((b - a) / N) / m
    h2 = 0.5 * h:  h6 = h / 6#
    u = a: v = Ya: i = 0
    Do
        T(i) = u:  y(i) = v
        i = i + 1:  If i > N Then Exit Sub
        For j = 1 To m
            P0 = EDO1_F(nf, u, v)
            u = u + h2
            Pm1 = EDO1_F(nf, u, v + h2 * P0)
            Pm2 = EDO1_F(nf, u, v + h2 * Pm1)
            u = u + h2
            p1 = EDO1_F(nf, u, v + h * Pm2)
            v = v + h6 * (P0 + 2# * (Pm1 + Pm2) + p1)
        Next j
    Loop
End Sub
Function EDO1_MERSON(ByVal nf As Integer, _
                     ByVal a As Double, _
                     ByVal b As Double, _
                     ByVal Ya As Double, _
                     ByVal N As Integer, _
                     ByVal m As Integer, _
                     ByRef T() As Double, _
                     ByRef y() As Double) As Double

    Dim i As Integer, j As Integer
    Dim h As Double, h2 As Double, h3 As Double
    Dim u As Double, v As Double
    Dim P0 As Double, p1 As Double, p2 As Double
    Dim P3 As Double, P4 As Double
    Dim E As Double, Emax As Double
    
    If N < 1 Or m < 1 Then Stop
    ReDim T(0 To N), y(0 To N)
    h = ((b - a) / N) / m
    h2 = 0.5 * h:  h3 = h / 3#
    u = a: v = Ya: i = 0
    Emax = 0: E = 0
    Do
        T(i) = u:  y(i) = v
        i = i + 1:  If i > N Then Exit Do
        For j = 1 To m
            P0 = h * EDO1_F(nf, u, v)
            p1 = h * EDO1_F(nf, u + h3, v + P0 / 3)
            p2 = 3 * h * EDO1_F(nf, u + h3, v + (P0 + p1) / 6)
            P3 = 4 * h * EDO1_F(nf, u + h2, v + (P0 + p2) / 8)
            P4 = h * EDO1_F(nf, u + h, v + (P0 - p2 + P3) / 2)
            u = u + h
            v = v + (P0 + P3 + P4) / 6
            E = (2 * P0 - 3 * p2 + 2 * P3 - P4) / 30
            If Abs(E) > Emax Then Emax = E
        Next j
    Loop
    EDO1_MERSON = Emax
End Function
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'      Système de deux Equations du premier ordre
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Sub EDO2_RK4(ByVal nf As Integer, _
             ByVal a As Double, _
             ByVal b As Double, _
             ByVal Ya As Double, _
             ByVal Za As Double, _
             ByVal N As Integer, _
             ByVal m As Integer, _
             ByRef T() As Double, _
             ByRef y() As Double, _
             ByRef z() As Double)

    Dim i As Integer, j As Integer
    Dim h As Double, h2 As Double, h6 As Double
    Dim u As Double, v As Double, w As Double
    Dim P0y As Double, Pm1y As Double, Pm2y As Double, P1y As Double
    Dim P0z As Double, Pm1z As Double, Pm2z As Double, P1z As Double
    
    If N < 1 Or m < 1 Then Stop
    ReDim T(0 To N), y(0 To N), z(0 To N)
    h = ((b - a) / N) / m
    h2 = 0.5 * h:    h6 = h / 6#
    u = a: v = Ya:     w = Za
    i = 0
    Do
        T(i) = u: y(i) = v: z(i) = w
        i = i + 1: If i > N Then Exit Sub
        For j = 1 To m
            Call EDO2_F(nf, u, v, w, P0y, P0z)
            u = u + h2
            Call EDO2_F(nf, u, v + h2 * P0y, w + h2 * P0z, Pm1y, Pm1z)
            Call EDO2_F(nf, u, v + h2 * Pm1y, w + h2 * Pm1z, Pm2y, Pm2z)
            u = u + h2
            Call EDO2_F(nf, u, v + h * Pm2y, w + h2 * Pm2z, P1y, P1z)
            v = v + h6 * (P0y + 2# * (Pm1y + Pm2y) + P1y)
            w = w + h6 * (P0z + 2# * (Pm1z + Pm2z) + P1z)
        Next j
    Loop
End Sub
' <><><><><><><><><><><><><><><><><><><><><><><><><>
'
'      Système de S Equations du premier ordre
'
' <><><><><><><><><><><><><><><><><><><><><><><><><>
Sub EDOS_RK4(ByVal nf As Integer, _
             ByVal a As Double, _
             ByVal b As Double, _
             ByVal N As Integer, _
             ByVal m As Integer, _
             ByVal s As Integer, _
             ByRef matY() As Double)

    Dim i As Integer, j As Integer, k As Integer
    Dim h As Double, h2 As Double, h6 As Double
    Dim u As Double, v() As Double, w() As Double
    Dim P0() As Double, Pm1() As Double, Pm2() As Double, p1() As Double
    
    If N < 1 Or m < 1 Or s < 3 Then Stop
    ReDim P0(1 To s), Pm1(1 To s), Pm2(1 To s), p1(1 To s)
    ReDim v(1 To s), w(1 To s)
    h = ((b - a) / N) / m
    h2 = 0.5 * h:  h6 = h / 6#
    u = matY(0, 0)
    For k = 1 To s
        v(k) = matY(0, k)
    Next k
    i = 0
    For i = 1 To N
        For j = 1 To m
            Call EDOS_F(nf, u, v(), P0())
            u = u + h2
            For k = 1 To s
                w(k) = v(k) + h2 * P0(k)
            Next k
            Call EDOS_F(nf, u, w(), Pm1())
            For k = 1 To s
                w(k) = v(k) + h2 * Pm1(k)
            Next k
            Call EDOS_F(nf, u, w(), Pm2())
            u = u + h2
            For k = 1 To s
                w(k) = v(k) + h * Pm2(k)
            Next k
            Call EDOS_F(nf, u, w(), p1())
            For k = 1 To s
                v(k) = v(k) + h6 * (P0(k) + 2# * (Pm1(k) + Pm2(k)) + p1(k))
            Next k
        Next j
        matY(0, i) = u
        For k = 1 To s
            matY(i, k) = v(k)
        Next k
    Next i
End Sub

