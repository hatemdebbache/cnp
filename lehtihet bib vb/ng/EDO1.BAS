Attribute VB_Name = "EDO1"
Option Explicit
Sub Euler(ByVal nf As Integer, _
          ByVal a As Double, _
          ByVal b As Double, _
          ByVal Ya As Double, _
          ByVal N As Integer, _
          ByRef t() As Double, _
          ByRef y() As Double, _
 Optional ByVal m As Integer = 1)
          
 Dim i As Integer, j As Integer
 Dim h As Double, u As Double, v As Double
 
 ' contrôle des données
 If N < 0 Or m < 1 Then Stop
 ' Début du programme
 ReDim t(N), y(N)
 h = ((b - a) / N) / m
 u = a
 v = Ya
 i = 0
 Do
    t(i) = u
    y(i) = v
    i = i + 1
    If i > N Then Exit Sub
    For j = 1 To m
        v = v + h * Fonc2(nf, u, v)
        u = u + h
    Next j
 Loop
End Sub
Sub Heun(ByVal nf As Integer, _
         ByVal a As Double, _
         ByVal b As Double, _
         ByVal Ya As Double, _
         ByVal N As Integer, _
         ByRef t() As Double, _
         ByRef y() As Double, _
Optional ByVal m As Integer = 1)
          
 Dim i As Integer, j As Integer
 Dim h As Double, u As Double, v As Double
 Dim h2 As Double, Po As Double
 
 ' contrôle des données
 If N < 0 Or m < 1 Then Stop
 ' Début du programme
 ReDim t(N), y(N)
 h = ((b - a) / N) / m
 h2 = 0.5 * h
 u = a
 v = Ya
 i = 0
 Do
    t(i) = u
    y(i) = v
    i = i + 1
    If i > N Then Exit Sub
    For j = 1 To m
        Po = Fonc2(nf, u, v)
        u = u + h
        v = v + h2 * (Po + Fonc2(nf, u, v + h * Po))
    Next j
 Loop
End Sub
Sub PointMilieu(ByVal nf As Integer, _
                ByVal a As Double, _
                ByVal b As Double, _
                ByVal Ya As Double, _
                ByVal N As Integer, _
                ByRef t() As Double, _
                ByRef y() As Double, _
       Optional ByVal m As Integer = 1)
          
 Dim i As Integer, j As Integer
 Dim h As Double, u As Double, v As Double
 Dim h2 As Double, Po As Double
 
 ' contrôle des données
 If N < 0 Or m < 1 Then Stop
 ' Début du programme
 ReDim t(N), y(N)
 h = ((b - a) / N) / m
 h2 = 0.5 * h
 u = a
 v = Ya
 i = 0
 Do
    t(i) = u
    y(i) = v
    i = i + 1
    If i > N Then Exit Sub
    For j = 1 To m
        Po = Fonc2(nf, u, v)
        v = v + h * Fonc2(nf, u + h2, v + h2 * Po)
        u = u + h
    Next j
 Loop
End Sub
Sub RK4(ByVal nf As Integer, _
        ByVal a As Double, _
        ByVal b As Double, _
        ByVal Ya As Double, _
        ByVal N As Integer, _
        ByRef t() As Double, _
        ByRef y() As Double, _
Optional ByVal m As Integer = 1)
          
 Dim i As Integer, j As Integer
 Dim h As Double, u As Double, v As Double
 Dim h2 As Double, Po As Double, P1 As Double
 Dim h6 As Double, Pm1 As Double, pm2 As Double
 
 ' contrôle des données
 If N < 0 Or m < 1 Then Stop
 ' Début du programme
 ReDim t(N), y(N)
 h = ((b - a) / N) / m
 h2 = 0.5 * h
 h6 = h / 6#
 u = a
 v = Ya
 i = 0
 Do
    t(i) = u
    y(i) = v
    i = i + 1
    If i > N Then Exit Sub
    For j = 1 To m
        Po = Fonc2(nf, u, v)
        u = u + h2
        Pm1 = Fonc2(nf, u, v + h2 * Po)
        pm2 = Fonc2(nf, u, v + h2 * Pm1)
        u = u + h2
        P1 = Fonc2(nf, u, v + h * pm2)
        v = v + h6 * (Po + P1 + 2 * Pm1 + 2 * pm2)
    Next j
 Loop
End Sub

