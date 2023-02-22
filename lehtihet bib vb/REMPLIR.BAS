Attribute VB_Name = "REMPLIR"
'
'  Calcul Numérique & Programmation
'  A/H.E. Lehtihet
'
'  Module   : << REMPLIR >>
'  Requis   :
'  Externes : Fonc(nf,x)
'
'  Descrip. : Remplissage de vecteurs/matrices
'
Option Explicit
Sub Vecteur1(ByVal a As Double, _
             ByVal b As Double, _
             ByVal k1 As Integer, _
             ByVal k2 As Integer, _
             ByRef x() As Double)

    Dim k As Integer, h As Double, v As Double
    
    If k2 <= k1 Then Stop
    h = (b - a) / (k2 - k1)
    v = a
    For k = k1 To k2
        x(k) = v
        v = v + h
    Next k
End Sub
Sub Vecteur2(ByVal nf As Integer, _
             ByVal a As Double, _
             ByVal b As Double, _
             ByVal k1 As Integer, _
             ByVal k2 As Integer, _
             ByRef x() As Double, _
             ByRef y() As Double)
 
    Dim k As Integer, h As Double, v As Double

    If k2 <= k1 Then Stop
    h = (b - a) / (k2 - k1)
    v = a
    For k = k1 To k2
        x(k) = v
        y(k) = Fonc(nf, v)
        v = v + h
    Next k
End Sub
Sub Vecteur1rnd(ByVal a As Double, _
                ByVal b As Double, _
                ByVal k1 As Integer, _
                ByVal k2 As Integer, _
                ByRef x() As Double)
 
    Dim k As Integer, d As Double

    If k2 < k1 Then Stop
    d = b - a
    For k = k1 To k2
        x(k) = a + d * Rnd
    Next k
End Sub
Sub Vecteur2rnd(ByVal nf As Integer, _
                ByVal a As Double, _
                ByVal b As Double, _
                ByVal k1 As Integer, _
                ByVal k2 As Integer, _
                ByRef x() As Double, _
                ByRef y() As Double)
   
    Dim k As Integer, d As Double, v As Double

    If k2 < k1 Then Stop
    d = b - a
    For k = k1 To k2
        v = a + d * Rnd
        x(k) = v
        y(k) = Fonc(nf, v)
    Next k
End Sub
Sub Vecteur2fix(ByVal nf As Integer, _
                ByVal k1 As Integer, _
                ByVal k2 As Integer, _
                ByRef x() As Double, _
                ByRef y() As Double)
   
    Dim k As Integer

    If k2 < k1 Then Stop
    For k = k1 To k2
        y(k) = Fonc(nf, x(k))
    Next k
End Sub
Sub MatriceBord(ByVal ik1 As Integer, _
                ByVal ik2 As Integer, _
                ByVal jk1 As Integer, _
                ByVal jk2 As Integer, _
                ByRef MAT() As Double, _
                ByVal vik1 As Double, _
                ByVal vik2 As Double, _
                ByVal vjk1 As Double, _
                ByVal vjk2 As Double, _
       Optional ByVal moy As Double = -1E+30)
 
    Dim i As Integer, j As Integer
    
    If ik2 <= ik1 Then Stop
    If jk2 <= jk1 Then Stop
    If moy <> -1E+30 Then
        If moy = 1E+30 Then moy = 0.25 * (vik1 + vik2 + vjk1 + vjk2)
        For i = ik1 To ik2
            For j = jk1 To jk2
                MAT(i, j) = moy
            Next j
        Next i
    End If
    For i = ik1 + 1 To ik2 - 1
        MAT(i, jk1) = vjk1
        MAT(i, jk2) = vjk2
    Next i
    For j = jk1 + 1 To jk2 - 1
        MAT(ik1, j) = vik1
        MAT(ik2, j) = vik2
    Next j
End Sub
